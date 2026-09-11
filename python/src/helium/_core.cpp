#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "HeI_Atom.h"
#include "HeI_Ric_Topbase.h"

double const_HI_A2s_1s = const_HI_A2s_1s_0;
double const_sigT = const_sigT_0;
double const_HeI_A2s_1s = const_HeI_A2s_1s_0;

using namespace std;

namespace py = pybind11;

namespace
{
    const int default_j = -9999;

    string normalize_data_path(string data_path)
    {
        if(data_path.empty()) data_path = "./Helium.vX/Helium.Data/";
        if(data_path[data_path.size()-1] != '/') data_path += "/";
        return data_path;
    }

    double frequency_to_ev(double nu)
    { return const_h*nu/const_e; }

    class HeliumAtomWrapper
    {
    public:
        HeliumAtomWrapper(int n_shells,
                          int j_resolved_shells,
                          int quadrupole_shells,
                          int intercombination_shells,
                          string data_path,
                          bool initialize_photoionization,
                          int message_level):
            gas_(),
            n_shells_(n_shells),
            j_resolved_shells_(j_resolved_shells),
            quadrupole_shells_(quadrupole_shells),
            intercombination_shells_(intercombination_shells),
            data_path_(normalize_data_path(data_path)),
            message_level_(message_level),
            photoionization_ready_(false),
            topbase_ready_(false)
        {
            if(n_shells < 2) throw invalid_argument("n_shells must be at least 2");
            if(j_resolved_shells < 0) throw invalid_argument("j_resolved_shells must be non-negative");
            if(quadrupole_shells < 0) throw invalid_argument("quadrupole_shells must be non-negative");
            if(intercombination_shells < 0) throw invalid_argument("intercombination_shells must be non-negative");

            gas_.init(n_shells,
                      j_resolved_shells,
                      quadrupole_shells,
                      intercombination_shells,
                      data_path_,
                      message_level);

            n_shells_ = gas_.Get_nShells();
            j_resolved_shells_ = gas_.Get_njres();
            quadrupole_shells_ = gas_.Get_n_added_Q_lines();
            intercombination_shells_ = gas_.Get_n_added_TS_lines();
            data_path_ = gas_.Get_data_path();

            if(initialize_photoionization) init_photoionization_rates();
        }

        int n_shells() const { return n_shells_; }
        int j_resolved_shells() const { return j_resolved_shells_; }
        int quadrupole_shells() const { return quadrupole_shells_; }
        int intercombination_shells() const { return intercombination_shells_; }
        string data_path() const { return data_path_; }
        int message_level() const { return message_level_; }
        int level_count() const { return gas_.Get_total_number_of_Levels(); }
        bool photoionization_ready() const { return photoionization_ready_; }

        py::dict summary() const
        {
            py::dict result;
            result["n_shells"] = n_shells_;
            result["j_resolved_shells"] = j_resolved_shells_;
            result["quadrupole_shells"] = quadrupole_shells_;
            result["intercombination_shells"] = intercombination_shells_;
            result["data_path"] = data_path_;
            result["level_count"] = gas_.Get_total_number_of_Levels();
            result["photoionization_ready"] = photoionization_ready_;
            return result;
        }

        py::dict level(int n, int l, int s, int j) const
        {
            const int resolved_j = resolve_j(n, l, s, j);
            return level_to_dict(n, l, s, resolved_j);
        }

        py::list levels() const
        {
            py::list result;

            for(int n=1; n<=n_shells_; n++)
            {
                for(int l=0; l<n; l++)
                    result.append(level_to_dict(n, l, 0, l));
            }

            for(int n=2; n<=n_shells_; n++)
            {
                if(n<=j_resolved_shells_)
                {
                    for(int l=0; l<n; l++)
                    {
                        if(l==0) result.append(level_to_dict(n, l, 1, 1));
                        else
                        {
                            for(int j=l-1; j<=l+1; j++)
                                result.append(level_to_dict(n, l, 1, j));
                        }
                    }
                }
                else
                {
                    for(int l=0; l<n; l++)
                        result.append(level_to_dict(n, l, 1, -10));
                }
            }

            return result;
        }

        py::dict transition(int upper_n,
                            int upper_l,
                            int upper_s,
                            int upper_j,
                            int lower_n,
                            int lower_l,
                            int lower_s,
                            int lower_j) const
        {
            const int uj = resolve_j(upper_n, upper_l, upper_s, upper_j);
            const int lj = resolve_j(lower_n, lower_l, lower_s, lower_j);

            if(lower_n > upper_n)
                throw invalid_argument("lower_n must not be larger than upper_n");

            const double a21 = gas_.Get_A(upper_n, upper_l, upper_s, uj,
                                          lower_n, lower_l, lower_s, lj);
            const double nu21 = gas_.Get_nu21(upper_n, upper_l, upper_s, uj,
                                              lower_n, lower_l, lower_s, lj);
            const double lambda21 = gas_.Get_lambda21(upper_n, upper_l, upper_s, uj,
                                                      lower_n, lower_l, lower_s, lj);

            py::dict result;
            result["upper_n"] = upper_n;
            result["upper_l"] = upper_l;
            result["upper_S"] = upper_s;
            result["upper_J"] = uj;
            result["lower_n"] = lower_n;
            result["lower_l"] = lower_l;
            result["lower_S"] = lower_s;
            result["lower_J"] = lj;
            result["a21_per_s"] = a21;
            result["frequency_hz"] = nu21;
            result["wavelength_cm"] = lambda21;
            result["wavelength_angstrom"] = lambda21*1.0e8;
            result["exists"] = (a21 > 0.0);
            return result;
        }

        double two_photon_rate(int upper_n,
                               int upper_l,
                               int upper_s,
                               int upper_j,
                               int lower_n,
                               int lower_l,
                               int lower_s,
                               int lower_j) const
        {
            const int uj = resolve_j(upper_n, upper_l, upper_s, upper_j);
            const int lj = resolve_j(lower_n, lower_l, lower_s, lower_j);

            if(lower_n!=1 || lower_l!=0 || lower_s!=0 || lj!=0)
                throw invalid_argument("two-photon rates are provided only to 1^1S_0");

            if(upper_n==2 && upper_l==0 && upper_s==0 && uj==0)
                return const_HeI_A2s_1s;
            if(upper_n==2 && upper_l==0 && upper_s==1 && uj==1)
                return const_HeI_A23s_1s;

            return 0.0;
        }

        void init_photoionization_rates()
        {
            gas_.init_photoionization_rates(message_level_);
            photoionization_ready_ = true;
            topbase_ready_ = true;
        }

        double photoionization_rate(int n, int l, int s, int j, double radiation_temperature)
        {
            validate_temperature(radiation_temperature, "radiation_temperature");
            ensure_photoionization_ready();
            const int resolved_j = resolve_j(n, l, s, j);
            return gas_.R_ic(n, l, s, resolved_j, radiation_temperature);
        }

        double recombination_rate(int n, int l, int s, int j, double radiation_temperature, double rho)
        {
            validate_temperature(radiation_temperature, "radiation_temperature");
            validate_temperature(rho, "rho");
            ensure_photoionization_ready();
            const int resolved_j = resolve_j(n, l, s, j);
            return gas_.R_ci(n, l, s, resolved_j, radiation_temperature, rho);
        }

        double detailed_balance_recombination_rate(int n, int l, int s, int j, double radiation_temperature)
        {
            validate_temperature(radiation_temperature, "radiation_temperature");
            ensure_photoionization_ready();
            const int resolved_j = resolve_j(n, l, s, j);
            return gas_.R_ci_DB(n, l, s, resolved_j, radiation_temperature);
        }

        double topbase_cross_section(int n, int l, int s, int j, double frequency_hz)
        {
            validate_frequency(frequency_hz);
            resolve_j(n, l, s, j);
            validate_topbase_level(n, l, s);
            load_topbase_data();
            return sig_ic_Topbase(n, l, s, frequency_hz);
        }

        py::array_t<double> topbase_cross_section_array(
            int n,
            int l,
            int s,
            int j,
            py::array_t<double, py::array::c_style | py::array::forcecast> frequencies_hz)
        {
            resolve_j(n, l, s, j);
            validate_topbase_level(n, l, s);
            load_topbase_data();

            return map_array(frequencies_hz, [&](double frequency) {
                validate_frequency(frequency);
                return sig_ic_Topbase(n, l, s, frequency);
            });
        }

        py::dict topbase_photoionization(
            int n,
            int l,
            int s,
            int j,
            py::array_t<double, py::array::c_style | py::array::forcecast> frequencies_hz)
        {
            const int resolved_j = resolve_j(n, l, s, j);
            const int index = gas_.Get_Level_index(n, l, s, resolved_j);
            const double threshold = gas_.Get_nu_ion(index);
            py::array_t<double> sigma = topbase_cross_section_array(n, l, s, resolved_j, frequencies_hz);

            py::dict result;
            result["n"] = n;
            result["l"] = l;
            result["S"] = s;
            result["J"] = resolved_j;
            result["frequency_hz"] = py::array(frequencies_hz);
            result["threshold_hz"] = threshold;
            result["sigma_cm2"] = sigma;
            return result;
        }

        void rescale(double alpha_scale, double electron_mass_scale)
        {
            validate_temperature(alpha_scale, "alpha_scale");
            validate_temperature(electron_mass_scale, "electron_mass_scale");
            gas_.rescale_gas(alpha_scale, electron_mass_scale);
        }

        void reset()
        { gas_.reset_gas(); }

    private:
        Gas_of_HeI_Atoms gas_;
        int n_shells_;
        int j_resolved_shells_;
        int quadrupole_shells_;
        int intercombination_shells_;
        string data_path_;
        int message_level_;
        bool photoionization_ready_;
        bool topbase_ready_;

        int resolve_j(int n, int l, int s, int j) const
        {
            validate_nls(n, l, s);

            if(s==0)
            {
                const int expected_j = l;
                if(j==default_j) return expected_j;
                if(j!=expected_j) throw invalid_argument("singlet levels require J = l");
                return j;
            }

            if(n<=j_resolved_shells_)
            {
                if(l==0)
                {
                    if(j==default_j) return 1;
                    if(j!=1) throw invalid_argument("triplet S levels require J = 1");
                    return j;
                }

                if(j==default_j) return l;
                if(j<l-1 || j>l+1)
                    throw invalid_argument("triplet levels require J = l-1, l, or l+1");
                return j;
            }

            if(j==default_j) return -10;
            if(j!=-10) throw invalid_argument("non-j-resolved triplet levels use J = -10");
            return j;
        }

        void validate_nls(int n, int l, int s) const
        {
            if(n < 1 || n > n_shells_)
                throw out_of_range("n is outside the initialized shell range");
            if(l < 0 || l >= n)
                throw out_of_range("l must satisfy 0 <= l < n");
            if(s!=0 && s!=1)
                throw invalid_argument("S must be 0 for singlet or 1 for triplet");
            if(s==1 && n<2)
                throw invalid_argument("triplet helium levels start at n = 2");
        }

        void validate_topbase_level(int n, int l, int s) const
        {
            bool supported = false;

            if(n<=10)
            {
                if(l==0) supported = true;
                else if(l==1)
                {
                    if(s==0 && n<=9) supported = true;
                    if(s==1 && n<=3) supported = true;
                }
                else if(l==2 && n<=9) supported = true;
            }

            if(!supported)
                throw invalid_argument("TOPbase cross sections are not available for this level");
        }

        void validate_temperature(double value, const string &name) const
        {
            if(!isfinite(value) || value <= 0.0)
                throw invalid_argument(name + " must be finite and positive");
        }

        void validate_frequency(double frequency_hz) const
        {
            if(!isfinite(frequency_hz) || frequency_hz <= 0.0)
                throw invalid_argument("frequency_hz must be finite and positive");
        }

        void ensure_photoionization_ready()
        {
            if(!photoionization_ready_) init_photoionization_rates();
        }

        void load_topbase_data()
        {
            if(topbase_ready_) return;
            load_all_Topbase_data(data_path_ + "TopBase_data/");
            topbase_ready_ = true;
        }

        py::dict level_to_dict(int n, int l, int s, int j) const
        {
            const int index = gas_.Get_Level_index(n, l, s, j);
            const int ground = gas_.Get_Level_index(1, 0, 0, 0);
            const double nu_ion = gas_.Get_nu_ion(index);
            const double nu_to_ground = (index == ground ? 0.0 : gas_.Get_nu_ul(index, ground));

            py::dict result;
            result["index"] = index;
            result["n"] = n;
            result["l"] = l;
            result["S"] = s;
            result["J"] = j;
            result["statistical_weight"] = gas_.Get_gw(index);
            result["ionization_frequency_hz"] = nu_ion;
            result["ionization_energy_ev"] = frequency_to_ev(nu_ion);
            result["frequency_to_ground_hz"] = nu_to_ground;
            result["energy_to_ground_ev"] = frequency_to_ev(nu_to_ground);
            return result;
        }

        template <typename Function>
        py::array_t<double> map_array(
            py::array_t<double, py::array::c_style | py::array::forcecast> input_array,
            Function function) const
        {
            py::buffer_info input_info = input_array.request();
            vector<py::ssize_t> shape(input_info.shape.begin(), input_info.shape.end());
            py::array_t<double> output(shape);
            py::buffer_info output_info = output.request();

            const double *input = static_cast<const double *>(input_info.ptr);
            double *values = static_cast<double *>(output_info.ptr);

            for(py::ssize_t i=0; i<input_info.size; i++)
                values[i] = function(input[i]);

            return output;
        }
    };

    HeliumAtomWrapper make_helium(int shells,
                                  int j_resolved_shells,
                                  int quadrupole_shells,
                                  int intercombination_shells,
                                  string data_path,
                                  bool initialize_photoionization,
                                  int message_level)
    {
        return HeliumAtomWrapper(shells,
                                 j_resolved_shells,
                                 quadrupole_shells,
                                 intercombination_shells,
                                 data_path,
                                 initialize_photoionization,
                                 message_level);
    }
}

PYBIND11_MODULE(_core, module)
{
    module.doc() = "Python bindings for the Helium-demo C++ atom routines.";
    module.attr("heI_two_photon_2s_1s_rate_per_s") = const_HeI_A2s_1s;
    module.attr("heI_triplet_two_photon_2s_1s_rate_per_s") = const_HeI_A23s_1s;

    py::class_<HeliumAtomWrapper>(module, "HeliumAtom")
        .def(py::init<int, int, int, int, string, bool, int>(),
             py::arg("n_shells") = 10,
             py::arg("j_resolved_shells") = 10,
             py::arg("quadrupole_shells") = 10,
             py::arg("intercombination_shells") = 10,
             py::arg("data_path") = "./Helium.vX/Helium.Data/",
             py::arg("initialize_photoionization") = true,
             py::arg("message_level") = -1)
        .def_property_readonly("n_shells", &HeliumAtomWrapper::n_shells)
        .def_property_readonly("j_resolved_shells", &HeliumAtomWrapper::j_resolved_shells)
        .def_property_readonly("quadrupole_shells", &HeliumAtomWrapper::quadrupole_shells)
        .def_property_readonly("intercombination_shells", &HeliumAtomWrapper::intercombination_shells)
        .def_property_readonly("data_path", &HeliumAtomWrapper::data_path)
        .def_property_readonly("message_level", &HeliumAtomWrapper::message_level)
        .def_property_readonly("level_count", &HeliumAtomWrapper::level_count)
        .def_property_readonly("photoionization_ready", &HeliumAtomWrapper::photoionization_ready)
        .def("summary", &HeliumAtomWrapper::summary)
        .def("level",
             &HeliumAtomWrapper::level,
             py::arg("n"),
             py::arg("l"),
             py::arg("S") = 0,
             py::arg("J") = default_j)
        .def("levels", &HeliumAtomWrapper::levels)
        .def("transition",
             &HeliumAtomWrapper::transition,
             py::arg("upper_n"),
             py::arg("upper_l"),
             py::arg("upper_S"),
             py::arg("upper_J"),
             py::arg("lower_n"),
             py::arg("lower_l"),
             py::arg("lower_S") = 0,
             py::arg("lower_J") = default_j)
        .def("two_photon_rate",
             &HeliumAtomWrapper::two_photon_rate,
             py::arg("upper_n"),
             py::arg("upper_l"),
             py::arg("upper_S"),
             py::arg("upper_J"),
             py::arg("lower_n") = 1,
             py::arg("lower_l") = 0,
             py::arg("lower_S") = 0,
             py::arg("lower_J") = default_j)
        .def("init_photoionization_rates", &HeliumAtomWrapper::init_photoionization_rates)
        .def("photoionization_rate",
             &HeliumAtomWrapper::photoionization_rate,
             py::arg("n"),
             py::arg("l"),
             py::arg("S"),
             py::arg("J"),
             py::arg("radiation_temperature"))
        .def("recombination_rate",
             &HeliumAtomWrapper::recombination_rate,
             py::arg("n"),
             py::arg("l"),
             py::arg("S"),
             py::arg("J"),
             py::arg("radiation_temperature"),
             py::arg("rho") = 1.0)
        .def("detailed_balance_recombination_rate",
             &HeliumAtomWrapper::detailed_balance_recombination_rate,
             py::arg("n"),
             py::arg("l"),
             py::arg("S"),
             py::arg("J"),
             py::arg("radiation_temperature"))
        .def("topbase_cross_section",
             py::overload_cast<int, int, int, int, double>(&HeliumAtomWrapper::topbase_cross_section),
             py::arg("n"),
             py::arg("l"),
             py::arg("S"),
             py::arg("J"),
             py::arg("frequency_hz"))
        .def("topbase_cross_section",
             py::overload_cast<int,
                               int,
                               int,
                               int,
                               py::array_t<double, py::array::c_style | py::array::forcecast>>(
                 &HeliumAtomWrapper::topbase_cross_section_array),
             py::arg("n"),
             py::arg("l"),
             py::arg("S"),
             py::arg("J"),
             py::arg("frequencies_hz"))
        .def("topbase_photoionization",
             &HeliumAtomWrapper::topbase_photoionization,
             py::arg("n"),
             py::arg("l"),
             py::arg("S"),
             py::arg("J"),
             py::arg("frequencies_hz"))
        .def("rescale",
             &HeliumAtomWrapper::rescale,
             py::arg("alpha_scale"),
             py::arg("electron_mass_scale"))
        .def("reset", &HeliumAtomWrapper::reset);

    module.def("helium",
               &make_helium,
               py::arg("shells") = 10,
               py::arg("j_resolved_shells") = 10,
               py::arg("quadrupole_shells") = 10,
               py::arg("intercombination_shells") = 10,
               py::arg("data_path") = "./Helium.vX/Helium.Data/",
               py::arg("initialize_photoionization") = true,
               py::arg("message_level") = -1);
}
