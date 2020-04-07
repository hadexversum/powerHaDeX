

# inputSeq1 = "KITEGKLVIWINGDKGYNGLAEVGKKFEKDTGIKVTVEHPDKLEEKFPQVAATGDGPDIIFWAHDRFGGYAQSGLLAEITPDKAFQDKLYPFTWDAVRYNGKLIAYPIAVEALSLIYNKDLLPNPPKTWEEIPALDKELKAKGKSALMFNLQEPYFTWPLIAADGGYAFKYENGKYDIKDVGVDNAGAKAGLTFLVDLIKNKHMNADTDYSIAEAAFNKGETAMTINGPWAWSNIDTSKVNYGVTVLPTFKGQPSKPFVGVLSAGINAASPNKELAKEFLENYLLTDEGLEAVNKDKPLGAVALKSYEEELAKDPRIAATMENAQKGEIMPNIPQMSAFWYAVRTAVINAASGRQTVDEALKDAQTRITK"
# inputSeq = strsplit(inputSeq1, "")[[1]]
# TempC = 15
# pH = 9
# molType = 'poly'
# pDread = 9
# ifCorr = 0
# kcHD <- fbmme_hd(inputSeq, pDread, TempC, molType, ifCorr)

fbmme_hd <- function(inputSeq, pDread, TempC, molType, ifCorr){


  Temp = TempC + 273.15; #convert to kelvin scale

  if(ifCorr == 1){
    pD = pDread+0.4
  }else{
    if(ifCorr == 0){
      pD = pDread
    }else{
      print('ifCorr must be 1 or 0!')
    }
  }

  D = 10^(-pD) #[D+]
  OD = 10^(pD - 15.05) #[OD-]

  R = 1.9858775 #gas constant(cal/(K*mol)); the value in spreadsheet is 1.987

  Ea_A = 14000 #unit: cal/mol
  Ea_B = 17000
  Ea_W = 19000

  Q7 = (1/Temp - 1/293)/R #1/(dT*R)
  Fta = exp(-Q7 * Ea_A)
  Ftb = exp(-Q7 * Ea_B)
  Ftw = exp(-Q7 * Ea_W)

  Ea_Asp = 1000 #unit: cal/mol
  Ea_Glu = 1083
  Ea_His = 7500

  pKc_Asp = -log10(10^(-4.48)*exp(-1*Ea_Asp*((1/Temp-1/278)/R)))
  pKc_Glu = -log10(10^(-4.93)*exp(-1*Ea_Glu*((1/Temp-1/278)/R)))
  pKc_His = -log10(10^(-7.42)*exp(-1*Ea_His*((1/Temp-1/278)/R)))

  Ka_poly = (10^(1.62))/60
  Ka_oligo = Ka_poly*2.34

  Kb_poly = (10^(10.05))/60
  Kb_oligo = Kb_poly*1.35

  Kw_poly = (10^(-1.5))/60
  Kw_oligo = Kw_poly*1.585


  if(molType == "poly"){
    Ka = Ka_poly
    Kb = Kb_poly
    Kw = Kw_poly
  }else
    if(molType =="oligo"){
      Ka = Ka_oligo
      Kb = Kb_oligo
      Kw = Kw_oligo
    }else{
      print('Unknown molType!(must be "poly" or "oligo")')
    }

  ###########################################################################
  ###reproduce FBMME_HD.XLS Table B:
  # # # ASP COOH D+	d
  # # # CYSTINE	C2	s
  # # # GLU COOH E+	e
  # # # PRO cis	Pc	p
  # # # N-Term	NT	n
  # # # C-Term	CT	c
  # # # -NHMe	NMe	m
  # # # H3CCO-	Ac	a

  Factors_list = strsplit('ARDdNCsGEeQHILKMFPpSTWYVncma', "")[[1]]
  Factors_values = matrix(c(0,0,0,0,
                          -0.590000000000000,-0.320000000000000,0.0767122542818456,0.220000000000000,
                          0,0,0,0, #Asp; calculate later
                          -0.900000000000000,-0.120000000000000,0.690000000000000,0.600000000000000,
                          -0.580000000000000,-0.130000000000000,0.490000000000000,0.320000000000000,
                          -0.540000000000000,-0.460000000000000,0.620000000000000,0.550000000000000,
                          -0.740000000000000,-0.580000000000000,0.550000000000000,0.460000000000000,
                          -0.220000000000000,0.218176047120386,0.267251569286023,0.170000000000000,
                          0,0,0,0, #Glu; calculate later
                          -0.600000000000000,-0.270000000000000,0.240000000000000,0.390000000000000,
                          -0.470000000000000,-0.270000000000000,0.0600000000000000,0.200000000000000,
                          0,0,0,0, #His; calculate later
                          -0.910000000000000,-0.590000000000000,-0.730000000000000,-0.230000000000000,
                          -0.570000000000000,-0.130000000000000,-0.576252727721677,-0.210000000000000,
                          -0.560000000000000,-0.290000000000000,-0.0400000000000000,0.120000000000000,
                          -0.640000000000000,-0.280000000000000,-0.00895484265292644,0.110000000000000,
                          -0.520000000000000,-0.430000000000000,-0.235859464059171,0.0631315866300978,
                          0,-0.194773472023435,0,-0.240000000000000,
                          0,-0.854416534276379,0,0.600000000000000,
                          -0.437992277698594,-0.388518934646472,0.370000000000000,0.299550285605933,
                          -0.790000000000000,-0.468073125742265,-0.0662579798400606,0.200000000000000,
                          -0.400000000000000,-0.440000000000000,-0.410000000000000,-0.110000000000000,
                          -0.410000000000000,-0.370000000000000,-0.270000000000000,0.0500000000000000,
                          -0.739022273362575,-0.300000000000000,-0.701934483299758,-0.140000000000000,
                          0,-1.32000000000000,0,1.62000000000000,
                          0,0,-1.80000000000000,0, #C-Term; col1 calculate later
                          0,0,0,0, #-NHMe; col1&3 calulate later
                          0,0.293000000000000,0,-0.197000000000000), ncol = 4, byrow=TRUE) #//above part same with that of FBMME_DH.XLS

  Factors_values[3,1] = log10(10^(-0.9-pD)/(10^(-pKc_Asp)+10^(-pD))+10^(0.9-pKc_Asp)/(10^(-pKc_Asp)+10^(-pD)))
  Factors_values[3,2] = log10(10^(-0.12-pD)/(10^(-pKc_Asp)+10^(-pD))+10^(0.58-pKc_Asp)/(10^(-pKc_Asp)+10^(-pD)))
  Factors_values[3,3] = log10(10^(0.69-pD)/(10^(-pKc_Asp)+10^(-pD))+10^(-0.3-pKc_Asp)/(10^(-pKc_Asp)+10^(-pD)))
  Factors_values[3,4] = log10(10^(0.6-pD)/(10^(-pKc_Asp)+10^(-pD))+10^(-0.18-pKc_Asp)/(10^(-pKc_Asp)+10^(-pD)))

  Factors_values[9,1] = log10(10^(-0.6-pD)/(10^(-pKc_Glu)+10^(-pD))+10^(-0.9-pKc_Glu)/(10^(-pKc_Glu)+10^(-pD)))
  Factors_values[9,2] = log10(10^(-0.27-pD)/(10^(-pKc_Glu)+10^(-pD))+10^(0.31-pKc_Glu)/(10^(-pKc_Glu)+10^(-pD)))
  Factors_values[9,3] = log10(10^(0.24-pD)/(10^(-pKc_Glu)+10^(-pD))+10^(-0.51-pKc_Glu)/(10^(-pKc_Glu)+10^(-pD)))
  Factors_values[9,4] = log10(10^(0.39-pD)/(10^(-pKc_Glu)+10^(-pD))+10^(-0.15-pKc_Glu)/(10^(-pKc_Glu)+10^(-pD)))

  Factors_values[12,1] = log10(10^(-0.8-pD)/(10^(-pKc_His)+10^(-pD))+10^(0-pKc_His)/(10^(-pKc_His)+10^(-pD)))
  Factors_values[12,2] = log10(10^(-0.51-pD)/(10^(-pKc_His)+10^(-pD))+10^(0-pKc_His)/(10^(-pKc_His)+10^(-pD)))
  Factors_values[12,3] = log10(10^(0.8-pD)/(10^(-pKc_His)+10^(-pD))+10^(-0.1-pKc_His)/(10^(-pKc_His)+10^(-pD)))
  Factors_values[12,4] = log10(10^(0.83-pD)/(10^(-pKc_His)+10^(-pD))+10^(0.14-pKc_His)/(10^(-pKc_His)+10^(-pD)))

  Factors_values[26,1] = log10(10^(0.05-pD)/(10^(-pKc_Glu)+10^(-pD))+10^(0.96-pKc_Glu)/(10^(-pKc_Glu)+10^(-pD)))

  Factors_values[27,1] = log10(135.5/(Ka*60))
  Factors_values[27,3] = log10(2970000000/(Kb*60))


  ###########################################################################
  ###calculate kc(H->D) of each position:

  Seq = c("n", inputSeq, "c")
  N = length(Seq) #residue numbers of the input peptide/protein plus 2(N- & C-term)
  kcHD = rep(0, N)

  for(i in 1:N){

    if(i==1 || Seq[i]=='P' || Seq[i]=='p' || Seq[i]=='a'){
      #do nothing, leave kc=0
    }else{
      ###Fa:
      if(Seq[i]=='P' || Seq[i]=='a' || i==1 || Seq[i]=='p' || i==N || i-1==1){
        Fa=0
      }else{
        j=which(Factors_list==Seq[i])
        k=which(Factors_list==Seq[i-1])

        if((i-2 <= 1 && i+1 == N && Seq[i-1] != 'a') && Seq[i] != 'm'){
          Fa=10^(Factors_values[j,1]+Factors_values[k,2]+Factors_values[25,2]+Factors_values[26,1])
        }else{
          if(i-2<=1 && Seq[i-1]!='a'){
            Fa=10^(Factors_values[j,1] + Factors_values[k,2] + Factors_values[25,2])
          }else{
            if( i+1 == N && Seq[i] != 'm'){
              Fa=10^(Factors_values[j,1] + Factors_values[k,2] + Factors_values[26,1])
            }else{
              Fa=10^(Factors_values[j,1] + Factors_values[k,2])
            }
          }
        }
      }
      ###Fb:
      if( Seq[i]=='P' || Seq[i]=='a' || i==1 || Seq[i]=='p' || i==N || i-1==1){
        Fb = 0
      }else{
        j = which(Factors_list == Seq[i])
        k = which(Factors_list == Seq[i-1])
        if((i - 2 <= 1 && i + 1 == N && Seq[i - 1] != 'a') && Seq[i] != 'm'){
          Fb = 10^(Factors_values[j, 3] + Factors_values[k, 4] + Factors_values[25, 4] + Factors_values[26, 3])
        }else{
          if(i - 2 <= 1 && Seq[i - 1] != 'a'){
            Fb=10^(Factors_values[j, 3] + Factors_values[k, 4] + Factors_values[25, 4])
          }else{
            if(i + 1 == N && Seq[i] != 'm'){
              Fb=10^(Factors_values[j, 3] + Factors_values[k, 4] + Factors_values[26, 3])
            }else{
              Fb=10^(Factors_values[j, 3] + Factors_values[k, 4])
            }
          }
        }
      }
    kcHD[i] = Fa*Ka*D*Fta + Fb*Kb*OD*Ftb + Fb*Kw*Ftw
    }
  }

  kcHD = kcHD[2:(N-1)]

  return(kcHD)
}

