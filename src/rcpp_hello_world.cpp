#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

const double gas_const = 1.9858775; // gas constant(cal/(K * mol)); the value in spreadsheet is 1.987
const double Ea_A = 14000; //unit: cal/mol
const double Ea_B = 17000;
const double Ea_W = 19000;
const double Ea_Asp = 960; //unit: cal/mol
const double Ea_Glu = 1083;
const double Ea_His = 7500;
const double Ka_poly = pow(10, 1.4) / 60.0;
const double Ka_oligo = Ka_poly * 2.34;

const double Kb_poly = pow(10, 9.87) / 60.0;
const double Kb_oligo = Kb_poly * 1.35;

const double Kw_poly = pow(10,-1.6) / 60.0;
const double Kw_oligo = Kw_poly * 1.585;

NumericMatrix get_constants_matrix(void) {
    NumericVector constants_dh = {0.0, 0.0, 0.0, 0.0, -0.59, -0.32, 0.0767122542818456, 0.22,
                                0.0, 0.0, 0.0, 0.0, // Asp; calculate later
                                -0.9, -0.12, 0.69, 0.6, -0.58, -0.13, 0.49, 0.32,
                                -0.54, -0.46, 0.62, 0.55, -0.74, -0.58, 0.55, 0.46,
                                -0.22, 0.218176047120386, 0.267251569286023, 0.17,
                                0.0, 0.0, 0.0, 0.0, // Glu; calculate later
                                -0.6, -0.27, 0.24, 0.39,
                                -0.47, -0.27, 0.06, 0.2,
                                0.0, 0.0, 0.0, 0.0, //His; calculate later
                                -0.91,-0.59,-0.73,-0.23,
                                -0.57,-0.13,-0.576252727721677,-0.21,
                                -0.56,-0.29,-0.040,0.12,
                                -0.64,-0.28,-0.00895484265292644,0.11,
                                -0.52,-0.43,-0.235859464059171,0.0631315866300978,
                                0.0,-0.194773472023435,0.0,-0.24,
                                0.0,-0.854416534276379,0.0,0.60,
                                -0.437992277698594,-0.388518934646472,0.37,0.299550285605933,
                                -0.79,-0.468073125742265,-0.0662579798400606,0.20,
                                -0.40,-0.44,-0.41,-0.11,
                                -0.41,-0.37,-0.27,0.050,
                                -0.739022273362575,-0.30,-0.701934483299758,-0.14,
                                0.0, -1.32, 0.0, 1.62,
                                0.0, 0.0, -1.8, 0.0, // C-Term; col1 calculate later
                                0.0, 0.0, 0.0, 0.0, // -NHMe; col1&3 calulate later
                                0.0, 0.293, 0.0, -0.197};
    NumericMatrix constants(28, 4);
    int j = 0;
    for (int i = 0; i < constants_dh.length(); ++i) {
        int col = i % 4;
        constants(j, col) = constants_dh[i];
        if(col == 0 && i != 0) {
            ++j;
        }
    }
    return constants;
}

// [[Rcpp::export]]
NumericMatrix get_constants(double pKc_Asp, double pKc_Glu, double pKc_His,
                            double Ka, double Kb, double pH)
{
    NumericMatrix constants = get_constants_matrix();
    constants(2,0) = log10(pow(10, -0.9-pH)/(pow(10, -pKc_Asp) + pow(10, -pH)) + pow(10, 0.9-pKc_Asp)/(pow(10, -pKc_Asp) + pow(10, -pH)));
    constants(2,1) = log10(pow(10, -0.12-pH)/(pow(10, -pKc_Asp) + pow(10, -pH)) + pow(10, 0.58-pKc_Asp)/(pow(10, -pKc_Asp) + pow(10, -pH)));
    constants(2,2) = log10(pow(10, 0.69-pH)/(pow(10, -pKc_Asp) + pow(10, -pH)) + pow(10, -0.3-pKc_Asp)/(pow(10, -pKc_Asp) + pow(10, -pH)));
    constants(2,3) = log10(pow(10, 0.6-pH)/(pow(10, -pKc_Asp) + pow(10, -pH)) + pow(10, -0.18-pKc_Asp)/(pow(10, -pKc_Asp) + pow(10, -pH)));

    constants(8,0) = log10(pow(10, -0.6-pH)/(pow(10, -pKc_Glu) + pow(10, -pH)) + pow(10, -0.9-pKc_Glu)/(pow(10, -pKc_Glu) + pow(10, -pH)));
    constants(8,1) = log10(pow(10, -0.27-pH)/(pow(10, -pKc_Glu) + pow(10, -pH)) + pow(10, 0.31-pKc_Glu)/(pow(10, -pKc_Glu) + pow(10, -pH)));
    constants(8,2) = log10(pow(10, 0.24-pH)/(pow(10, -pKc_Glu) + pow(10, -pH)) + pow(10, -0.51-pKc_Glu)/(pow(10, -pKc_Glu) + pow(10, -pH)));
    constants(8,3) = log10(pow(10, 0.39-pH)/(pow(10, -pKc_Glu) + pow(10, -pH)) + pow(10, -0.15-pKc_Glu)/(pow(10, -pKc_Glu) + pow(10, -pH)));

    constants(11,0) = log10(pow(10, -0.8-pH)/(pow(10, -pKc_His) + pow(10, -pH)) + pow(10, 0.0-pKc_His)/(pow(10, -pKc_His) + pow(10, -pH)));
    constants(11,1) = log10(pow(10, -0.51-pH)/(pow(10, -pKc_His) + pow(10, -pH)) + pow(10, 0.0-pKc_His)/(pow(10, -pKc_His) + pow(10, -pH)));
    constants(11,2) = log10(pow(10, 0.8-pH)/(pow(10, -pKc_His) + pow(10, -pH)) + pow(10, -0.1-pKc_His)/(pow(10, -pKc_His) + pow(10, -pH)));
    constants(11,3) = log10(pow(10, 0.83-pH)/(pow(10, -pKc_His) + pow(10, -pH)) + pow(10, 0.14-pKc_His)/(pow(10, -pKc_His) + pow(10, -pH)));

    constants(25,0) = log10(pow(10, 0.05-pH)/(pow(10, -pKc_Glu) + pow(10, -pH)) + pow(10, 0.96-pKc_Glu)/(pow(10, -pKc_Glu) + pow(10, -pH)));

    constants(26,0) = log10(135.5/(Ka * 60));
    constants(26,2) = log10(2970000000.0/(Kb * 60));

    return constants;
}

// [[Rcpp::export]]
CharacterVector fbmme_dh(CharacterVector sequence, double ph = 9.0, double temperature = 15.0,
                         String mol_type = "poly")
{
    const double kelvin_temp = temperature + 273.15;
    const double D = pow(10, -ph); // [D+]
    const double OD = pow(10, ph-14.17); // [OD-]
    const double Q7 = (1 / kelvin_temp-1 / 293) / gas_const; // 1/(dT * R)
    const double Fta = exp(-Q7 * Ea_A);
    const double Ftb = exp(-Q7 * Ea_B);
    const double Ftw = exp(-Q7 * Ea_W);

    const double pKc_Asp = -log10(pow(10, -3.87) * exp(-1 * Ea_Asp * ((1.0/kelvin_temp-1.0/278.0) / gas_const)));
    const double pKc_Glu = -log10(pow(10, -4.33) * exp(-1 * Ea_Glu * ((1.0/kelvin_temp-1.0/278.0) / gas_const)));
    const double pKc_His = -log10(pow(10, -7.0) * exp(-1 * Ea_His * ((1.0/kelvin_temp-1.0/278.0) / gas_const)));

    const double Ka = (mol_type == String("poly")) ? Ka_poly : Ka_oligo;
    const double Kb = (mol_type == String("poly")) ? Kb_poly : Kb_oligo;
    const double Kw = (mol_type == String("poly")) ? Kw_poly : Kw_oligo;

    // create the numeric matrix here
    sequence.push_front(String("n"));
    sequence.push_back(String("c"));

    int N = sequence.length();
    NumericVector kcDH (N);

    CharacterVector amino_acids {"A", "R", "D", "d", "N", "C", "s", "G", "E", "e",
                                 "Q", "H", "I", "L", "K", "M", "F", "P", "p", "S",
                                 "T", "W", "Y", "V", "n", "c", "m", "a"};
    // for(int i = 0; i < N; ++i) {
    //     if(i == 0 || sequence[i] == String("P") || sequence[i] == String("p") || sequence[i] == String("a")) {
    //         continue;
    //     } else {
    //         // Fa
    //         String current = which(Factors_list == Seq[i])
    //         String previous = which(Factors_list == Seq[i-1])
    //
    //
    //     }
    // }
    kcDH.erase(N-1);
    kcDH.erase(0);
    return amino_acids;
}


// //Fa:
// if(Seq[i] == 'P' || Seq[i] == 'a' || i == 1 || Seq[i] == 'p' || i == N || i-1 == 1){
//     Fa = 0
//
// }else{
//
//     if(i-2 <= 1 && i + 1 == N && Seq[i-1] != 'a' && Seq[i] != 'm'){
//         Fa = 10^(Factors_values[j,1] + Factors_values[k,2] + Factors_values[25,2] + Factors_values[26,1])
//
//     }else{
//         if(i-2 <= 1 && Seq[i-1] != 'a'){
//             Fa = 10^(Factors_values[j,1] + Factors_values[k,2] + Factors_values[25,2])
//
//         }else{
//             if(i + 1 == N && Seq[i] != 'm'){
//                 Fa = 10^(Factors_values[j,1] + Factors_values[k,2] + Factors_values[26,1])
//
//             }else{
//                 Fa = 10^(Factors_values[j,1] + Factors_values[k,2])
//             }
//         }
//     }
// }
//
// //Fb:
// if(Seq[i] == 'P' || Seq[i] == 'a' || i == 1 || Seq[i] == 'p' || i == N || i-1 == 1){
//     Fb = 0
//
// }else{
//     j = which(Factors_list == Seq[i])
//     k = which(Factors_list == Seq[i-1])
//
//     if(i-2 <= 1 && i + 1 == N && Seq[i-1] != 'a' && Seq[i] != 'm'){
//         Fb = 10^(Factors_values[j,3] + Factors_values[k,4] + Factors_values[25,4] + Factors_values[26,3])
//
//     }else{
//         if(i-2 <= 1 && Seq[i-1] != 'a'){
//             Fb = 10^(Factors_values[j,3] + Factors_values[k,4] + Factors_values[25,4])
//
//         }else{
//             if (i + 1 == N && Seq[i] != 'm'){
//                 Fb = 10^(Factors_values[j,3] + Factors_values[k,4] + Factors_values[26,3])
//             }else{
//                 Fb = 10^(Factors_values[j,3] + Factors_values[k,4])
//             }
//         }
//     }
//     kcDH[i] = Fa * Ka * D * Fta + Fb * Kb * OD * Ftb + Fb * Kw * Ftw
// }
// }
