//============================================================================================
//
// Simple demo for use of the Helium atom class. The Helium setup was developed for detailed
// cosmological recombination calculations, where singlet and triplet states, intercombination
// lines, quadrupole transitions, photoionization, and recombination rates all enter.
//
// References:
// [1] https://ui.adsabs.harvard.edu/abs/2011MNRAS.412..748C/abstract
// [2] https://ui.adsabs.harvard.edu/abs/2008A%26A...485..377R/abstract
// [3] https://ui.adsabs.harvard.edu/abs/2007ApJS..170..251D/abstract
//
//============================================================================================
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "HeI_Atom.h"
#include "HeI_Ric_Topbase.h"

double const_HI_A2s_1s = const_HI_A2s_1s_0;
double const_sigT = const_sigT_0;
double const_HeI_A2s_1s = const_HeI_A2s_1s_0;

using namespace std;

namespace
{
    struct LevelID
    {
        int n;
        int l;
        int s;
        int j;
        string label;
    };

    double frequency_to_ev(double nu)
    { return const_h*nu/const_e; }

    int level_index(const Gas_of_HeI_Atoms &helium, const LevelID &level)
    { return helium.Get_Level_index(level.n, level.l, level.s, level.j); }

    void print_atom_summary(const Gas_of_HeI_Atoms &helium)
    {
        cout << "Neutral helium example" << '\n';
        cout << "  shells=" << helium.Get_nShells()
             << "  j-resolved up to n=" << helium.Get_njres()
             << "  levels=" << helium.Get_total_number_of_Levels()
             << '\n';
        cout << "  data_path=" << helium.Get_data_path() << '\n';
        cout << "  added intercombination lines="
             << (helium.are_add_TS_lines_loaded() ? "yes" : "no")
             << "  added quadrupole lines="
             << (helium.are_add_Q_lines_loaded() ? "yes" : "no")
             << '\n';
    }

    void print_level_summary(const Gas_of_HeI_Atoms &helium, const LevelID &level)
    {
        const int index = level_index(helium, level);
        const int ground = helium.Get_Level_index(1, 0, 0, 0);
        const double nu_ion = helium.Get_nu_ion(index);
        const double nu_to_ground = (index == ground ? 0.0 : helium.Get_nu_ul(index, ground));

        cout << "  " << level.label
             << "  (n,l,S,J)=(" << level.n << "," << level.l << ","
             << level.s << "," << level.j << ")"
             << "  gw=" << helium.Get_gw(index)
             << "  E_ion=" << frequency_to_ev(nu_ion) << " eV"
             << "  DeltaE_1s2=" << frequency_to_ev(nu_to_ground) << " eV"
             << '\n';
    }

    void print_transition(const Gas_of_HeI_Atoms &helium,
                          const LevelID &upper,
                          const LevelID &lower,
                          const string &label)
    {
        const double A21 = helium.Get_A(upper.n, upper.l, upper.s, upper.j,
                                        lower.n, lower.l, lower.s, lower.j);
        const double nu21 = helium.Get_nu21(upper.n, upper.l, upper.s, upper.j,
                                            lower.n, lower.l, lower.s, lower.j);
        const double lambda21 = helium.Get_lambda21(upper.n, upper.l, upper.s, upper.j,
                                                    lower.n, lower.l, lower.s, lower.j);

        cout << "  " << label
             << ": A21=" << A21 << " 1/s"
             << "  nu21=" << nu21 << " Hz"
             << "  lambda21=" << lambda21 << " cm"
                 << '\n';
    }

    void print_two_photon_rate(const string &label, double rate)
    {
        cout << "  " << label
             << ": A2gamma=" << rate << " 1/s"
             << '\n';
    }

    void print_rate_summary(Gas_of_HeI_Atoms &helium, const LevelID &level, double Tg)
    {
        const int index = level_index(helium, level);

        cout << "  " << level.label
             << " at Tg=" << Tg << " K"
             << ": R_ic=" << helium.R_ic(index, Tg) << " 1/s"
             << "  R_ci=" << helium.R_ci(index, Tg) << " cm^3/s"
             << "  R_ci_DB=" << helium.R_ci_DB(index, Tg) << " cm^3/s"
             << '\n';
    }

    void print_topbase_cross_sections(const Gas_of_HeI_Atoms &helium,
                                      const LevelID &level,
                                      const vector<double> &threshold_factors)
    {
        const int index = level_index(helium, level);
        const double nu_ion = helium.Get_nu_ion(index);

        cout << "  " << level.label
             << "  threshold=" << nu_ion << " Hz"
             << '\n';

        for(size_t i=0; i<threshold_factors.size(); i++)
        {
            const double factor = threshold_factors[i];
            const double nu = factor*nu_ion;

            cout << "    nu/nu_ion=" << factor
                 << "  sigma_TopBase=" << sig_ic_Topbase(level.n, level.l, level.s, nu)
                 << " cm^2"
                 << '\n';
        }
    }
}

//============================================================================================
int main()
{
    const string helium_data_path = "./Helium.vX/Helium.Data/";
    const int shells = 10;
    const int j_resolved_shells = 10;
    const int quadrupole_shells = 10;
    const int intercombination_shells = 10;
    const int message_level = -1;

    Gas_of_HeI_Atoms helium(shells,
                            j_resolved_shells,
                            quadrupole_shells,
                            intercombination_shells,
                            helium_data_path,
                            message_level);

    helium.init_photoionization_rates(message_level);

    const LevelID ground = {1, 0, 0, 0, "1^1S_0"};
    const LevelID singlet_2s = {2, 0, 0, 0, "2^1S_0"};
    const LevelID singlet_2p = {2, 1, 0, 1, "2^1P_1"};
    const LevelID singlet_3d = {3, 2, 0, 2, "3^1D_2"};
    const LevelID triplet_2s = {2, 0, 1, 1, "2^3S_1"};
    const LevelID triplet_2p = {2, 1, 1, 1, "2^3P_1"};

    cout << scientific << setprecision(6);
    cout << "Helium atom demo using Helium.vX" << '\n';
    print_atom_summary(helium);

    cout << "\nSelected levels" << '\n';
    print_level_summary(helium, ground);
    print_level_summary(helium, singlet_2s);
    print_level_summary(helium, singlet_2p);
    print_level_summary(helium, triplet_2s);
    print_level_summary(helium, triplet_2p);

    cout << "\nSelected bound-bound and forbidden transitions" << '\n';
    print_transition(helium, singlet_2p, ground, "He I resonance line 2^1P_1 -> 1^1S_0");
    print_transition(helium, triplet_2p, ground, "He I intercombination 2^3P_1 -> 1^1S_0");
    print_transition(helium, singlet_3d, ground, "He I quadrupole 3^1D_2 -> 1^1S_0");

    cout << "\nSelected two-photon channels" << '\n';
    print_two_photon_rate("He I 2^1S_0 -> 1^1S_0", const_HeI_A2s_1s);
    print_two_photon_rate("He I 2^3S_1 -> 1^1S_0", const_HeI_A23s_1s);

    cout << "\nPhotoionization and recombination rates" << '\n';
    print_rate_summary(helium, singlet_2s, 3000.0);
    print_rate_summary(helium, singlet_2p, 3000.0);
    print_rate_summary(helium, triplet_2s, 3000.0);

    cout << "\nTOPbase photoionization cross sections" << '\n';
    print_topbase_cross_sections(helium, ground, {1.001, 1.1, 2.0, 5.0});
    print_topbase_cross_sections(helium, singlet_2p, {1.001, 1.5, 3.0});
    print_topbase_cross_sections(helium, triplet_2s, {1.001, 1.5, 3.0});

    const double nu0 = helium.Get_nu21(singlet_2p.n, singlet_2p.l, singlet_2p.s, singlet_2p.j,
                                       ground.n, ground.l, ground.s, ground.j);
    helium.rescale_gas(1.01, 1.0);
    const double nu_scaled = helium.Get_nu21(singlet_2p.n, singlet_2p.l, singlet_2p.s, singlet_2p.j,
                                             ground.n, ground.l, ground.s, ground.j);

    cout << "\nRescaling example" << '\n';
    cout << "  alpha_scale=1.01, me_scale=1.0"
         << "  2^1P_1 -> 1^1S_0 frequency ratio=" << nu_scaled/nu0
         << "  original=" << nu0 << " Hz"
         << "  scaled=" << nu_scaled << " Hz"
         << '\n';

    helium.reset_gas();

    return 0;
}

//============================================================================================
//============================================================================================
