// Sperry Model -- C++ version for GNU compilers

#include <stdio.h>
#include <cmath> // math utility functions
#include <string> // the C++ String Class, easier to deal with than char arrays for this application
#include <ctime> // timers for performance testing

#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>

#define FIO_PRECISION 12

#define STAGE_ID_NONE 0 // use this mode to skip all of the stage code and just run based on param sheet settings (a "normal" run, the default config)
#define STAGE_ID_HIST_OPT 1
#define STAGE_ID_HIST_STRESS 2
#define STAGE_ID_FUT_OPT 3
#define STAGE_ID_FUT_STRESS 4
#define STAGE_ID_FUT_STRESS_NOACCLIM 5

#define MAX_SUMMARY_COLS 121
#define MAX_SUMMARY_ROWS 2001

std::string stageNames[STAGE_ID_FUT_STRESS_NOACCLIM + 1];

const long bufferYears = 5; // throw out this many years when calculating optimal BA:GA to avoid the field capacity starting condition influencing the results
double dummyDouble = 0.0;

#define PARAMFILE_MAXROWS 2001
#define PARAMFILE_MAXCOLS 101
#define DATAFILE_MAXROWS 2000001
#define DATAFILE_MAXCOLS 101

std::string paramCells[PARAMFILE_MAXROWS][PARAMFILE_MAXCOLS];
std::string nameTable[PARAMFILE_MAXROWS][PARAMFILE_MAXCOLS]; // needs to be same dimensions as the parameter file "cells"
double dataCells[DATAFILE_MAXROWS][DATAFILE_MAXCOLS]; // up to 2000k rows and 100 columns of input data
double finalOutCells[MAX_SUMMARY_ROWS][MAX_SUMMARY_COLS];

long GSCells[101][11]; // contains the growing season start/end days
                       // a much simpler class for the full C++ version which just wraps access to the dataSet[][] array in the dsheet.Cells command
                       // exists to preserve the code legacy of the model, which is written in VBA for Excel
class CPP_IO_Handler
{
public:
   CPP_IO_Handler() {};
   double& Cells(long row, long col);
   double& fCells(long row, long col);
};

// returns a reference to the "cell" in question for reading or writing
// for IO code equivalency with original VBA version
double& CPP_IO_Handler::Cells(long row, long col)
{
   double outDbl = 0.0;
   if (row >= 0 && row < DATAFILE_MAXROWS)
   {
      if (col >= 0 && col < DATAFILE_MAXCOLS)
      {
         return dataCells[row][col];
      }
   }
   return dummyDouble;
}

double& CPP_IO_Handler::fCells(long row, long col)
{
   double outDbl = 0.0;
   if (row >= 0 && row < MAX_SUMMARY_ROWS)
   {
      if (col >= 0 && col < MAX_SUMMARY_COLS)
      {
         return finalOutCells[row][col];
      }
   }
   return dummyDouble;
}

CPP_IO_Handler dSheet;

class CSVRow
{
public:
   std::string const& operator[](std::size_t index) const
   {
      return m_data[index];
   }
   std::size_t size() const
   {
      return m_data.size();
   }
   void readNextRow(std::istream& str)
   {
      std::string         line;
      std::getline(str, line);

      std::stringstream   lineStream(line);
      std::string         cell;

      m_data.clear();
      while (std::getline(lineStream, cell, ','))
      {
         m_data.push_back(cell);
      }
      // This checks for a trailing comma with no data after it.
      if (!lineStream && cell.empty())
      {
         // If there was a trailing comma then add an empty element.
         m_data.push_back("");
      }
   }
private:
   std::vector<std::string>    m_data;
};

std::istream& operator>>(std::istream& str, CSVRow& data)
{
   data.readNextRow(str);
   return str;
}

double fAbs(double in)
{
   if (in < 0) return 0 - in;
   else
      return in;
}

double getValueFromNameDbl(std::string name)
{
   // search the names array for this string
   for (int x = 0; x < PARAMFILE_MAXROWS; x++)
   {
      for (int y = 0; y < PARAMFILE_MAXCOLS; y++)
      {
         if (nameTable[x][y] == name)
            return atof(paramCells[x][y].c_str());
      }
   }
   return 0.0;
}

std::string getValueFromNameStr(std::string name)
{
   // search the names array for this string
   for (int x = 0; x < PARAMFILE_MAXROWS; x++)
   {
      for (int y = 0; y < PARAMFILE_MAXCOLS; y++)
      {
         if (nameTable[x][y] == name)
            return paramCells[x][y];
      }
   }
   return "";
}

long getValueFromNameLng(std::string name)
{
   // search the names array for this string
   for (int x = 0; x < PARAMFILE_MAXROWS; x++)
   {
      for (int y = 0; y < PARAMFILE_MAXCOLS; y++)
      {
         if (nameTable[x][y] == name)
            return atoi(paramCells[x][y].c_str());
      }
   }
   return 0;
}

long getRowFromName(std::string name)
{
   // search the names array for this string
   for (int x = 0; x < PARAMFILE_MAXROWS; x++)
   {
      for (int y = 0; y < PARAMFILE_MAXCOLS; y++)
      {
         if (nameTable[x][y] == name)
            return x;
      }
   }
   return 0;
}

long getColFromName(std::string name)
{
   // search the names array for this string
   for (int x = 0; x < PARAMFILE_MAXROWS; x++)
   {
      for (int y = 0; y < PARAMFILE_MAXCOLS; y++)
      {
         if (nameTable[x][y] == name)
            return y;
      }
   }
   return 0;
}

class ModelProgram
{
public:
   double *p_kr, *p_erh, *p_er, *p_krh, *p_tkr, *p_ter;
   double *p_elayer, *p_prhizo;
   double *p_jmatrix;
   double *p_kroot;
   double e, br, cr, n[6], a[6], x, ksatrh[6], kmaxrh[6], ksatr[6], kmaxr[6], ksh;
   double bs, cs, bl, cl, ksats, ksatl, el[100001], es[100001], pcrits, pcritl, ps, tel[100001], tes[100001];
   double pinc, einc, pd[6], ksatp, rsatp, ksatroot, vp;
   double s, del, eps, epsx, sthresh, sum, sums, olds, kmin;
   double kr[6][100001], erh[6][100001], er[6][100001], krh[6][100001], tkr[6][100001], ter[6][100001];
   double elow, ehigh, plow, prh[6], estart, efinish, flow, dfrdpr, frt;
   double klow, khigh, kupper, klower, dfrhdprh[7], dfrdprh[7], dfrhdpr[7], pr;
   double pcritrh[6], pcritr[6], threshold, kfactor;
   double vv[7], aamax, jmatrix[7][7], dum, indx[7], func[7];
   double elayer[6][100001], prhizo[6][100001], pstem[100001], proot[100001], kleaf[100001], kminleaf;
   double kstem[100001], kminstem, kroot[6][100001], kminroot[6], eplant[100001], md, ecritsystem, pcritsystem; /*phigh*/
   long phigh;

   // [HNT] store virgin curves
   double el_v[100001], es_v[100001], er_v[6][100001], kr_v[6][100001];
   // [/HNT]

   std::string failspot, layerfailure[6], tlayerfailure[6], setting, refilling, ground, soilred;
   double dp, kplant[100001], kminplant, toplayer;
   double dedp[100001], pleaf[100001], dedpf[1000001], p1, p2, pl;
   double plold, predawn, dedpl, prtarget;
   double kmaxs, pstarget, kmaxl, pltarget, rhizsat, rootsat, sum2, beta, depth[6], rhizok;
   double patm, vpd, dedplzero;
   double gcmd, dpmax, gmax, vpdsat;
   long pstop, halt, failure, k, check, t, j, it, f, dd, tmax, tnm, layer[6], tlayer[6], test, layers, d, z, p, pmax, unknowns, i, imax, ii, ll, total;
   double kplantold, frac, eincdef, aspect;
   double length[6], radius[6], vol, vertdistance[11], depthmax, shallow, coef;
   double rhizotarg, rplant, rstem, rootr, rleaf, rhizor, rrhizofrac, kinc, kplantmax;
   double vgterm;
   double ci, ca, gcanw[100001], gcanc[100001], marker, var, psyn[100001], transpiration;
   double dedpmax, amaxfrac[100001];
   double dpa[100001], lavpdmd;
   double psynact, psynmax, dpamax;
   double lavpd[100001], grad, gha, numerator, denominator, lambda, emiss, airtemp, leaftemp[100001], rabs;
   double laperba, par, qmax, vmax25, kc25, ko25, comp25, theta, wind, leafwidth;
   double comp, vmax, kc, ko, je, jc, eplantl[1000001], lavpdc1, lavpdsum, lavpdh, gmaxl, cin[1000001];
   double cinc, jmax, jmax25, havmax, hdvmax, svvmax, hajmax, hdjmax, svjmax, jact, pleafv[1000001], maxvpd;
   double kloss[1000001], klossv[1000001], maxkloss;
   long chalk, skip, psynstop, totalv;
   double lsc, timestep, leafpercent, dbh, height, pgrav, ffc, pground, grounddistance;
   double water[6], fc[6], baperga, layerflow, soilredist[6], groundwater, deficit, rain, sumrain, swclimit[6], store, pend;
   double kkmax[6], thetasat[6], groundflow, rday25, rday[100001], waterold, waternew, waterchange, soilf[6];
   double drainage, gwflow, thetac;
   double lat, longitude, slope, slopeasp, lai, xang;
   double fet, et, sm, lc, tsn, sindec, dec, cosdec, tim, tod, coszen, zen, cosaz, az;
   double m, sp, sb, sd, st, cloud, obssolar, fcd, kbe, kbezero, mleafang;
   double rad, k1, t1, told, t2, kd, qd, qds, qdt, qb, qbt, qsc, qsh, qsl, laisl, laish, parsh, parsl, parbottom;
   double nirsh, nirsl, sshade, ssun, sbottom, ssunb, ssund, sref, ppfd, ea, eac, la, lg;
   long jd, idiot, o;
   double mdsh, cincsh, lavpdmdsh, gcmdsh, psynactsh, transpirationsh;
   long haltsh;
   double leaftempsh[100001], lavpdsh[100001], rdaysh[100001];
   double gcanwsh[100001], gcancsh[100001], cinsh[100001], psynmaxsh, psynsh[100001], amaxfracsh[100001];
   double transpirationtree, anetsh, anettree, atree, anet, rabssoil, alt, tsncorr;
   double rha, rhs, soilevap, soilabssol, rough, zdispl, xh, zh, us, mdensair, soilep, emission, soiltemp;
   std::string night, sevap, pet, rainsim;
   bool nightBool;
   double tau, minwind;
   double thetafc[6], thetafracres[6], thetafracfc[6];
   long runmean, weird, sign, xx, runs, ticks;
   double rmean, cutoff, dedplmin, amaxmax, lightcurv;
   double emd, leaftmd, leaftshmd, lavpdshmd, gcanwmd, gcancmd, rdaymd, psynmd[100001], psynmaxmd;
   double cinmd, gcanwshmd, gcancshmd, rdayshmd, psynshmd[100001], cinshmd, psynmaxshmd;
   double prinitial, dpasun, lightcomp, fieldcapfrac, rockfrac;
   double initialthreshold, inp[101];
   double kref, lscref, pdref, leafpercentref, kpday1, kxday1, runoff;
   long kmaxits, bst, startbst, scenario;
   bool kmaxset, leafpercentset;
   long year_cur, year_start; //not to be confused with the year array index, which is year_cur - year_start
   long yearVal; // the temporary variable where we hold the year read from the sheet

   double dpamin;

   long iter_Counter; //how many iterations of this data set have we run? For finding the supply curve
   bool iter_useAreaTable;
   bool iter_yearsAsCount; //instead of trying to keep track of year offsets (year index = start - cur),
                           //use years as dataset ids directly (year index = cur year)

   bool iter_gwEnable;
   double iter_gwInc; //how much should be increment the ground water every iteration? Smaller = higher resolution supply curve
   double iter_gwStart; //Some data sets may not respond to increasing ground water until a high threshold, so can start with a minimum level and increment from there
   double iter_gwEnd;
   bool iter_runSupplyCurve;
   double iter_gwDist; //we//ll keep track of the gw dist separately when running iterations, then override the grounddistance variable with this value
   long iter_ddOutMod; //output offset so we can keep multiple iterations of data

   bool iter_ffcEnable;
   double iter_ffc;
   double iter_ffcStart;
   double iter_ffcEnd;
   double iter_ffcInc;
   double iter_gwRunning;

   bool iter_bagaEnable;
   double iter_baga;
   double iter_bagaStart;
   double iter_bagaEnd;
   double iter_bagaInc;
   double iter_bagaRef;
   double iter_bagaCutoff;

   // these are for the new BAGA optimization method... reset at start of each stage, NOT on years or iterations
   double iter_refK;
   double iter_prevMeanK;
   double iter_prevMinK;
   double iter_prevBAGA;
   long iter_code; // used to override the normal iteration behavior if we want to walk back a half step

   double iter_prevUnder_BAGA;
   double iter_prevUnder_MeanK;
   double iter_prevUnder_minK;
   long iter_prevUnder_code; // 0 = 10% mean PLC, 1 = 85% max PLC, 2 = LAI

   double iter_prevOver_BAGA;
   double iter_prevOver_MeanK;
   double iter_prevOver_minK;
   long iter_prevOver_code; // 1 = 10% mean PLC, 2 = 85% max PLC, 3 = LAI

   double iter_override_next_baga; // set the iter_code to override and then put in the next BAGA, to totally ignore the normal routine
   double iter_override_max_baga; // calculate the maximum possible BA:GA based on the LAI and store it

   bool iter_bisect_overHalf; // is the value over BAmax/2? if YES we need to subtract 1 from the iteration counter, because we did 1 extra initial iteration at the same size range
   long iter_bisect_underCount;
   long iter_bisect_exponent; // store the exponent so we can increment
   bool iter_bisect_foundOver; // have we found an "over" value? if not we'll keep iterating by 1/8 of the max until we do
                               // end

   std::string raining;
   bool rainEnabled;

   double treeToPhotoLAI;

   bool useGSData;


   long rowD, colD, rowLR, colLR;

   long totalDataPoints, lastDataRow, lastDataCol;
   double peakMatchTSNCorr, peakMatchTau;

   long dColYear, dColDay, dColTime, dColSolar, dColWind, dColRain, dColTAir, dColTSoil, dColD;
   //Data column positions - NOTE THAT THEY ARE OFFSET FROM colD, the starting data column (to avoid wasting io array space)
   long dColF_p1, dColF_p2, dColF_p3, dColF_p4, dColF_p5, dColF_predawn, dColF_P, dColF_E, dColF_Gw, dColF_laVPD, dColF_leaftemp, dColF_ANet,
   dColF_s1m2, dColF_ci, dColF_PPFD, dColF_S_P, dColF_S_E, dColF_S_Gw, dColF_S_laVPD, dColF_S_leaftemp,
   dColF_S_Anet, dColF_S_s1m2, dColF_S_ci, dColF_S_PPFD;

   long dColF_T_E, dColF_T_ANet, dColF_T_s1m2,
   dColF_T_pcrit, dColF_T_Ecrit, dColF_CP_Pstem, dColF_CP_Proot, dColF_CP_kstem, dColF_CP_kleaf, dColF_CP_kplant,
   dColF_CP_kxylem, dColF_CP_kroot1, dColF_CP_kroot2, dColF_CP_kroot3, dColF_CP_kroot4, dColF_CP_kroot5, dColF_CP_krootAll,
   dColF_CP_Eroot1, dColF_CP_Eroot2, dColF_CP_Eroot3, dColF_CP_Eroot4, dColF_CP_Eroot5, dColF_CP_Empty1, dColF_CP_Empty2,
   dColF_End_watercontent, dColF_End_waterchange, dColF_End_rain, dColF_End_gwater, dColF_End_E, dColF_End_drainage,
   dColF_End_soilEvap, dColF_End_ET, dColF_End_ANet, dColF_End_input, dColF_End_PLCplant, dColF_End_PLCxylem,
   dColF_End_runoff;

   long dColF_GS_year, dColF_GS_input, dColF_GS_Anet, dColF_GS_E, dColF_GS_PLCp, dColF_GS_PLCx, dColF_GS_kPlant, dColF_GS_kXylem, dColF_GS_ET;

   long gs_ar_years[100];
   long gs_ar_starts[100];
   long gs_ar_ends[100];
   long growSeasonCount;

   double gs_ar_input[100]; // all of these are arrays of size 100 - but this should never matter as long as the gs_yearIndex never exceeds 99
   double gs_ar_Anet[100];
   double gs_ar_E[100];
   double gs_ar_PLCp[100];
   double gs_ar_PLCx[100];
   double gs_ar_kPlant[100];
   double gs_ar_kXylem[100];
   double gs_ar_ET[100];
   double gs_ar_PLC85[100];
   double gs_ar_PLCSum[100];
   double gs_ar_PLCSum_N[100];
   double gs_ar_waterInitial[100];
   double gs_ar_waterFinal[100];

   double gs_ar_waterInitial_GS[100];
   double gs_ar_waterFinal_GS[100];
   double gs_ar_waterInput_GS[100];

   double gs_ar_waterInitial_OFF[100];
   double gs_ar_waterFinal_OFF[100];
   double gs_ar_waterInput_OFF[100];

   long gs_ar_nrFailConverge[100];
   double gs_ar_nrFailConverge_Water[100];
   double gs_ar_nrFailConverge_WaterMax[100];
   long gs_ar_nrFailThreshold[100];

   double gs_ar_kPlantMean[100];
   long gs_ar_kPlantMean_N[100];

   double gs_ar_cica[100];
   long gs_ar_cica_N[100];

   double gs_ar_Aci[100];
   double gs_ar_AnetDay[100];

   bool isNewYear;
   long gs_yearIndex; //this is a counter from 0 (for the first year) indicating how many years have passed
                      //get the actual year from gs_ar_years(gs_yearIndex)
                      //this avoids having to make year an input column in the model -- will just count how many years have passed when running
   long gs_prevDay;
   bool gs_inGrowSeason;
   bool gs_doneFirstDay; //done all first day of grow season calculations?

   long yearCounter;

   long stage_ID; // see below for details on stage IDs
                  // these never get erased, so they can be stored between run stages
                  //double stage_OptBAGA = 0.0; // maybe just use 1 var for this?
   double stage_OptHistBAGA = 0.0;
   double stage_OptFutBAGA = 0.0;
   double stage_KmaxFut = 0.0;
   double stage_LeafResistPerFut = 0.0;
   double stage_CO2Fut = 0.0;
   double stage_Hist_refK = 0.0;
   double stage_Fut_refK = 0.0;

   // the minimum info necessary to identify any historical run
   std::string name_region = "";
   std::string name_site = "";
   std::string name_species = "";
   // and future
   std::string name_scen = "";
   std::string name_model = "";

   const double pi = 3.14159;
   const double sbc = 0.0000000567; //'stefan boltzman constant in W m-2 K-4
   const double sha = 29.3; //'specific heat of air in J mol-1C-1
   const double gas = 8.3144598; //'universal gas constant J mol-1K-1
   const double oa = 0.21; //'mole fraction of o2
   const double solar = 1362; //'solar constant W m-2
                              //'Const tau = 0.65 'clear sky transmissivity, CN p. 173
   const double absolar = 0.5; //'absorptivity of solar for leaves
   const double abspar = 0.8; //'absorptivity of par for leaves
   const double absnir = 0.2; //'absorptivity of near infrared for leaves

   

   double rvg(double &x) //'gives soil Y in MPa from soil theta/thetasat=x
   {
      //rvg = (x ^ (1 / (1 - 1 / n(z))) + 1) ^ (1 / n(z)) / (x ^ (1 / (n(z) - 1)) * a(z))
      double aa = pow((pow(x, (1 / (1 - 1 / n[z]))) + 1), (1 / n[z]));
      double bb = (pow(x, (1 / (n[z] - 1))) * a[z]);
      return aa / bb;
   }

   double swc(double &x) //'gives soil water content, theta/thetasat from soil water potential in MPa=x
   {
      //swc = (1 / (1 + (a(z) * x) ^ n(z))) ^ (1 - 1 / n(z))
      return pow((1 / (1 + pow((a[z] * x), n[z]))), (1 - 1 / n[z]));
   }


   double vg(double &x) //'the van genuchten function for soil k
   {
      //vp = 1 / ((a(z) * x) ^ n(z) + 1)
      //vg = kmaxrh(z) * vp ^ ((n(z) - 1) / (2 * n(z))) * ((1 - vp) ^ ((n(z) - 1) / n(z)) - 1) ^ 2
      vp = 1 / (pow((a[z] * x), n[z]) + 1);
      return kmaxrh[z] * pow(vp, ((n[z] - 1) / (2 * n[z]))) * pow((pow((1 - vp), ((n[z] - 1) / n[z])) - 1), 2);
   }

   void trapzdvg(double &p1, double &p2, double &s, long &t) //integrates van genuchten
   {
      if (t == 1)
      {
         s = 0.5 * (p2 - p1) * (vg(p1) + vg(p2));
         it = 1;
      }
      else
      {
         tnm = it;
         del = (p2 - p1) / tnm;
         x = p1 + 0.5 * del;
         sum = 0;
         for (j = 1; j <= it; j++)
         {
            sum = sum + vg(x);
            x = x + del;
         }
         s = 0.5 * (s + (p2 - p1) * sum / tnm);
         it = 2 * it;
      }
   }

   void qtrapvg(double &p1, double &p2, double &s) //evaluates accuracy of van genuchten integration
   {
      eps = 0.001; //fractional refinement threshold for integral
      tmax = 70; //limit of calls to trapzdvg
      olds = -1; //starting point unlikely to satisfy if statement below
      for (t = 1; t <= tmax; t++)
      {
         trapzdvg(p1, p2, s, t);
         if (std::abs(s - olds) < eps * std::abs(olds))
            return;
         olds = s;
      }
   }

   //'ROOT BLOCK OF SUBROUTINES

   double wbr(double &x) //the weibull function for root elements
   {
      return ksatr[z] * exp(-(pow((x / br), cr)));
      //wbr = ksatr(z) * Exp(-((x / br) ^ cr))
   }

   void trapzdwbr(double &p1, double &p2, double &s, long &t) //integrates root element z weibull
   {
      if (t == 1)
      {
         s = 0.5 * (p2 - p1) * (wbr(p1) + wbr(p2));
         it = 1;
      }
      else
      {
         tnm = it;
         del = (p2 - p1) / tnm;
         x = p1 + 0.5 * del;
         sum = 0;
         for (j = 1; j <= it; j++)
         {
            sum = sum + wbr(x);
            x = x + del;
         }
         s = 0.5 * (s + (p2 - p1) * sum / tnm);
         it = 2 * it;
      }
      //[HNT]
      //testcounterIntRoot = testcounterIntRoot + 1
      //[\HNT]
   }

   void qtrapwbr(double &p1, double &p2, double &s) //'evaluates accuracy of root element z integration
   {
      olds = -1; //'starting point unlikely to satisfy if statement below
      for (t = 1; t <= f; t++)
      {
         trapzdwbr(p1, p2, s, t);
         if (std::abs(s - olds) <= epsx * std::abs(olds))
            return;
         olds = s;
      }
   }

   //'STEM BLOCK OF SUBROUTINES

   double wbs(double &x) //the weibull function for stem element
   {
      //wbs = ksats * Exp(-((x / bs) ^ cs))
      return ksats * exp(-(pow((x / bs), cs)));
   }

   void trapzdwbs(double &p1, double &p2, double &s, long &t) //integrates root element z weibull
   {
      if (t == 1)
      {
         s = 0.5 * (p2 - p1) * (wbs(p1) + wbs(p2));
         it = 1;
      }
      else
      {
         tnm = it;
         del = (p2 - p1) / tnm;
         x = p1 + 0.5 * del;
         sum = 0;
         for (j = 1; j <= it; j++)
         {
            sum = sum + wbs(x);
            x = x + del;
         }
         s = 0.5 * (s + (p2 - p1) * sum / tnm);
         it = 2 * it;
      }
   }

   void qtrapwbs(double &p1, double &p2, double &s) //evaluates accuracy of root element z integration
   {
      olds = -1; //starting point unlikely to satisfy if statement below
      for (t = 1; t <= f; t++)
      {
         trapzdwbs(p1, p2, s, t);
         if (std::abs(s - olds) <= epsx * std::abs(olds))
            return;
         olds = s;
      }
   }

   //'LEAF BLOCK OF SUBROUTINES

   double wbl(double &x) //the weibull function for leaf element
   {
      //wbl = ksatl * Exp(-((x / bl) ^ cl));
      return ksatl * exp(-(pow((x / bl), cl)));
   }

   void trapzdwbl(double &p1, double &p2, double &s, long &t) //integrates root element z weibull
   {
      if (t == 1)
      {
         s = 0.5 * (p2 - p1) * (wbl(p1) + wbl(p2));
         it = 1;
      }
      else
      {
         tnm = it;
         del = (p2 - p1) / tnm;
         x = p1 + 0.5 * del;
         sum = 0;
         for (j = 1; j <= it; j++)
         {
            sum = sum + wbl(x);
            x = x + del;
         }
         s = 0.5 * (s + (p2 - p1) * sum / tnm);
         it = 2 * it;
      }
   }

   void qtrapwbl(double &p1, double &p2, double &s) //evaluates accuracy of root element z integration
   {
      olds = -1; //starting point unlikely to satisfy if statement below
      for (t = 1; t <= f; t++)
      {
         trapzdwbl(p1, p2, s, t);
         if (std::abs(s - olds) <= epsx * std::abs(olds))
            return;
         olds = s;
      }
   }

   bool locateRanges()
   {
      // set the parameter and nametable filenames
      std::string paramFileName = "parameters.csv";
      std::string nameTableFileName = "nametable.csv";

      std::cout << "Reading parameters from " << paramFileName << std::endl;

      std::ifstream paramFile(paramFileName);
      if (!paramFile.is_open())
      {
         std::cout << "UNRECOVERABLE: Failed to open parameter file " << paramFileName << std::endl;
         return false;
      }

      // load the string values of all the parameters into an array
      // in most cases these need to be converted to doubles when called on:
      // double x = atof(stringVar.c_str());
      long rowCount = -1;
      CSVRow row;
      while (paramFile >> row)
      {
         rowCount = rowCount + 1;
         if (rowCount >= PARAMFILE_MAXROWS)
            break;

         for (int rC = 0; rC < row.size(); rC++)
         {
            if (rowCount < PARAMFILE_MAXROWS && rowCount >= 0 && rC < PARAMFILE_MAXCOLS && rC >= 0)
            {
               paramCells[rowCount + 1][rC + 1] = row[rC];
               //if (row[rC] != "" && row[rC] != "\r")
               //std::cout << "Reading parameters. " << row[rC] << std::endl;
            }
         }
         //std::cout << "4th Element(" << row[3] << ")\n";
      }

      // load the names of all the cells from the parameter sheet into a another array
      // this allows us to maintain the "name-based" access used in the readin function
      // this is a bit of a hack -- in theory it allows us to avoid referencing cell coords here,
      // but it also means that the "nameTable" needs to be maintained if anything is added/removed
      // can always access the param sheets directly with paramCells[row][col]
      std::cout << "Reading name table from " << nameTableFileName << std::endl;

      std::ifstream namesFile(nameTableFileName);
      if (!namesFile.is_open())
      {
         std::cout << "UNRECOVERABLE: Failed to open nametable file " << nameTableFileName << std::endl;
         return false;
      }

      rowCount = -1;
      while (namesFile >> row)
      {
         rowCount = rowCount + 1;
         if (rowCount >= PARAMFILE_MAXROWS)
            break;

         for (int rC = 0; rC < row.size(); rC++)
         {
            if (rowCount < PARAMFILE_MAXROWS && rowCount >= 0 && rC < PARAMFILE_MAXCOLS && rC >= 0)
            {
               nameTable[rowCount + 1][rC + 1] = row[rC];
               //if (row[rC] != "" && row[rC] != "\r")
               //std::cout << "Reading  names. " << row[rC] << std::endl;
            }
         }
         //std::cout << "4th Element(" << row[3] << ")\n";
      }

      rowD = getRowFromName("datamarker_hourly") + 1;
      colD = getColFromName("datamarker_hourly") - 1;

      rowLR = getRowFromName("inputmarker_layer"); // not a typo -- the header marker for data is one above where it should be, but the header for layers is not .. the conventions here are a mess
      colLR = getColFromName("inputmarker_layer") - 1;

      std::cout << "Name Table Test, rowLR = " << rowLR << std::endl;

      // an ugly hack to setup the column order
      dColYear = 1;
      dColDay = dColYear + 1;
      dColTime = dColDay + 1;
      dColSolar = dColTime + 1;
      dColRain = dColSolar + 1;
      dColWind = dColRain + 1;
      dColTAir = dColWind + 1;
      dColTSoil = dColTAir + 1;
      dColD = dColTSoil + 1;
      dColF_p1 = dColD + 1;
      dColF_p2 = dColF_p1 + 1;
      dColF_p3 = dColF_p2 + 1;
      dColF_p4 = dColF_p3 + 1;
      dColF_p5 = dColF_p4 + 1;
      dColF_predawn = dColF_p5 + 1;
      dColF_P = dColF_predawn + 1;
      dColF_E = dColF_P + 1;
      dColF_Gw = dColF_E + 1;
      dColF_laVPD = dColF_Gw + 1;
      dColF_leaftemp = dColF_laVPD + 1;
      dColF_ANet = dColF_leaftemp + 1;
      dColF_s1m2 = dColF_ANet + 1;
      dColF_ci = dColF_s1m2 + 1;
      dColF_PPFD = dColF_ci + 1;
      dColF_S_P = dColF_PPFD + 1;
      dColF_S_E = dColF_S_P + 1;
      dColF_S_Gw = dColF_S_E + 1;
      dColF_S_laVPD = dColF_S_Gw + 1;
      dColF_S_leaftemp = dColF_S_laVPD + 1;
      dColF_S_Anet = dColF_S_leaftemp + 1;
      dColF_S_s1m2 = dColF_S_Anet + 1;
      dColF_S_ci = dColF_S_s1m2 + 1;
      dColF_S_PPFD = dColF_S_ci + 1;
      dColF_T_E = dColF_S_PPFD + 1;
      dColF_T_ANet = dColF_T_E + 1;
      dColF_T_s1m2 = dColF_T_ANet + 1;
      dColF_T_pcrit = dColF_T_s1m2 + 1;
      dColF_T_Ecrit = dColF_T_pcrit + 1;
      dColF_CP_Pstem = dColF_T_Ecrit + 1;
      dColF_CP_Proot = dColF_CP_Pstem + 1;
      dColF_CP_kstem = dColF_CP_Proot + 1;
      dColF_CP_kleaf = dColF_CP_kstem + 1;
      dColF_CP_kplant = dColF_CP_kleaf + 1;
      dColF_CP_kxylem = dColF_CP_kplant + 1;
      dColF_CP_kroot1 = dColF_CP_kxylem + 1;
      dColF_CP_kroot2 = dColF_CP_kroot1 + 1;
      dColF_CP_kroot3 = dColF_CP_kroot2 + 1;
      dColF_CP_kroot4 = dColF_CP_kroot3 + 1;
      dColF_CP_kroot5 = dColF_CP_kroot4 + 1;
      dColF_CP_krootAll = dColF_CP_kroot5 + 1;
      dColF_CP_Eroot1 = dColF_CP_krootAll + 1;
      dColF_CP_Eroot2 = dColF_CP_Eroot1 + 1;
      dColF_CP_Eroot3 = dColF_CP_Eroot2 + 1;
      dColF_CP_Eroot4 = dColF_CP_Eroot3 + 1;
      dColF_CP_Eroot5 = dColF_CP_Eroot4 + 1;
      dColF_CP_Empty1 = dColF_CP_Eroot5 + 1;
      dColF_CP_Empty2 = dColF_CP_Empty1 + 1;
      dColF_End_watercontent = dColF_CP_Empty2 + 1;
      dColF_End_waterchange = dColF_End_watercontent + 1;
      dColF_End_rain = dColF_End_waterchange + 1;
      dColF_End_gwater = dColF_End_rain + 1;
      dColF_End_E = dColF_End_gwater + 1;
      dColF_End_drainage = dColF_End_E + 1;
      dColF_End_soilEvap = dColF_End_drainage + 1;
      dColF_End_ET = dColF_End_soilEvap + 1;
      dColF_End_ANet = dColF_End_ET + 1;
      dColF_End_input = dColF_End_ANet + 1;
      dColF_End_PLCplant = dColF_End_input + 1;
      dColF_End_PLCxylem = dColF_End_PLCplant + 1;
      dColF_End_runoff = dColF_End_PLCxylem + 1;

      dColF_GS_year = 1;//dColF_End_runoff + 1;
      dColF_GS_input = dColF_GS_year + 1;
      dColF_GS_Anet = dColF_GS_input + 1;
      dColF_GS_E = dColF_GS_Anet + 1;
      dColF_GS_PLCp = dColF_GS_E + 1;
      dColF_GS_PLCx = dColF_GS_PLCp + 1;
      dColF_GS_kPlant = dColF_GS_PLCx + 1;
      dColF_GS_kXylem = dColF_GS_kPlant + 1;
      dColF_GS_ET = dColF_GS_kXylem + 1;

      // TODO HACK get the region, site, spec names now -- these aren't in the current nameTable files so we'll do it manually
      // fix this to use the nameTable when we have time to re-run all the acclimations... or write a script to update all the nametables
      // note: this is specific to the BA/GA optimization
      name_region = paramCells[18][2];
      name_site = paramCells[19][2];
      name_species = paramCells[17][2];

      name_scen = paramCells[21][2];
      name_model = paramCells[20][2];
      std::cout << "Read site and scenario names (region | site | species | scenario | model): " << name_region << " | " << name_site << " | " << name_species << " | " << name_scen << " | " << name_model << std::endl;

      return true;
   }

   CSVRow dataHeaderRow;
   CSVRow summaryHeaderRow;

   const long maxYears = 90;

   void readGSSheet()
   {
      memset(GSCells, 0, sizeof(GSCells));

      std::string GSFileName = "";

      if (stage_ID == STAGE_ID_NONE)
      {
         GSFileName = "seasonlimits.csv";
      }

      if (useGSData && GSFileName != "") // only if we actually selected a file, and we're using the GS data files in the first place
      {

         //std::cout << "Reading growing season limits." << std::endl;

         std::ifstream dataFile(GSFileName);
         if (!dataFile.is_open())
            std::cout << "FAILED to open GS RANGE file " << GSFileName << std::endl;

         std::cout << "Opened GS RANGE file " << GSFileName << std::endl;

         long yearCount = 0; // how many unique years have we encountered? 
         long curYear = 0;
         long readYear = 0;

         int rowCount = -1;
         CSVRow row;
         while (dataFile >> row)
         {
            rowCount = rowCount + 1;
            if (rowCount > 100)
               break;

            if (rowCount + 1 == 1)
            {
               // header row
               //headerRow = row;
               ; // do nothing, we don't need the header row
                 // note that the data starts on row 2 as a result
            }
            else if (rowCount + 1 > 1)
            {
               readYear = std::atol(row[0].c_str());
               if (readYear > curYear)
               {
                  curYear = readYear;
                  yearCount++;
                  if (yearCount > maxYears)
                  {
                     std::cout << "Read GS range for " << yearCount - 1 << " years, quitting read." << std::endl;
                     break;
                  }
                  std::cout << "Reading GS range year " << curYear << " (" << yearCount << "/" << maxYears << ")" << std::endl;
               }

               for (int rC = 0; rC < row.size(); rC++)
               {
                  if (rowCount < 100 && rowCount >= 0 && rC < 10 && rC >= 0)
                     GSCells[rowCount + 1][rC + 1] = std::atol(row[rC].c_str()); // load array as double
                                                                                 // all data i/o is in double
               }
            }
            //std::cout << "4th Element(" << row[3] << ")\n";
         }
         std::cout << "Finished GS RANGE read." << std::endl;
      }
      else
      {
         std::cout << "Skipping GS RANGE load -- Years will be treated as independent." << std::endl;
      }
   }

   void readDataSheet()
   {
      // try to find a header file
      bool foundHeaderFile = false;
      std::string headerFileName = "dataheader.csv";

      //std::cout << "Reading data header." << std::endl;

      std::ifstream headerFile(headerFileName);
      if (!headerFile.is_open())
         std::cout << "FAILED to open DATA HEADER file " << headerFileName << std::endl;
      else
      {
         std::cout << "Reading DATA HEADER file " << headerFileName << std::endl;
         if (headerFile >> dataHeaderRow) // should be the only row in the file
         {
            foundHeaderFile = true;
            //std::cout << headerRow << std::endl;
         }
      }

      bool foundSumHeaderFile = false;
      headerFileName = "sumheader.csv";

      //std::cout << "Reading summary header." << std::endl;

      std::ifstream sumHeaderFile(headerFileName);
      if (!sumHeaderFile.is_open())
         std::cout << "FAILED to open SUMMARY HEADER file " << headerFileName << std::endl;
      else
      {
         std::cout << "Reading SUMMARY HEADER file " << headerFileName << std::endl;
         if (sumHeaderFile >> summaryHeaderRow) // should be the only row in the file
         {
            foundSumHeaderFile = true;
            //std::cout << headerRow << std::endl;
         }
      }

      // otherwise, we'll load whatever head exists in the weather file later

      memset(dataCells, 0, sizeof(dataCells));

      std::string dataFileName = "dataset.csv";

      if (stage_ID == STAGE_ID_NONE)
      {
         dataFileName = "dataset.csv";
      }

      //std::cout << "Reading data set." << std::endl;

      std::ifstream dataFile(dataFileName);
      if (!dataFile.is_open())
         std::cout << "FAILED to open DATA file " << dataFileName << std::endl;

      std::cout << "Reading DATA file " << dataFileName << std::endl;

      // only want to read 30 years of data
      // long maxYears = 30; // change this to increase limit // moved to global
      // for now we'll take the first 30 we encounter
      long yearCount = 0; // how many unique years have we encountered? 
      long curYear = 0;
      long readYear = 0;

      int rowCount = -1;
      CSVRow row;
      while (dataFile >> row)
      {
         rowCount = rowCount + 1;
         if (rowCount >= DATAFILE_MAXROWS) // ~maximum size of an Excel spreadsheet
            break;

         if (rowCount + 2 == rowD && !foundHeaderFile)
         {
            // header row
            dataHeaderRow = row;
         }
         else if (rowCount + 2 > rowD)
         {
            // check the years
            readYear = std::atol(row[0].c_str());
            if (readYear > curYear)
            {
               curYear = readYear;
               yearCount++;
               if (yearCount > maxYears)
               {
                  std::cout << "Read in " << yearCount - 1 << " years, quitting read." << std::endl;
                  break;
               }
               std::cout << "Reading year " << curYear << " (" << yearCount << "/" << maxYears << ")" << std::endl;
            }

            for (int rC = 0; rC < row.size(); rC++)
            {
               if (rowCount < DATAFILE_MAXROWS && rowCount >= 0 && rC < DATAFILE_MAXCOLS && rC >= 0)
                  dataCells[rowCount + 1][rC + 1] = atof(row[rC].c_str()); // load array as double
                                                                           // all data i/o is in double
               else
                  long breakpoint = 1137;
            }
         }
         //std::cout << "4th Element(" << row[3] << ")\n";
      }
      std::cout << "Finished DATA read." << std::endl;
   }

   void readin() //'inputs and calculates all parameters at the start
   {
      std::cout << "INIT: Setting up model parameters. (year count = " << gs_yearIndex << ")" << std::endl;

      treeToPhotoLAI = getValueFromNameDbl("i_treeToPhotoLAI"); // these names come from the "nametable" and correspond to cells in the parameters csv

      lightcomp = getValueFromNameDbl("i_lightComp"); //Cells(14, 21) //light compensation point in ppfd
      runmean = 1; //Cells(6, 9) //running mean for profit maximization
      cutoff = 1.1; //Cells(7, 9) //cutoff for stopping dpamax search
      minwind = 0.4515; //m s-1'minimum wind threshold
      f = 70; //itmax for trapzd routine used for xylem only
      epsx = 0.0001; //fractional acceptable error for trapzd routine in xylem
      sthresh = 1.0001; //acceptable error for e integral for xylem

      br = getValueFromNameDbl("i_br"); //getInputValue(wbRange, 3, 1) //Cells(4, 2) // = br //weibull b for each root element
      cr = getValueFromNameDbl("i_cr"); //getInputValue(wbRange, 4, 1) //Cells(5, 2) // = cr //weibull c for each root element
      bs = getValueFromNameDbl("i_bs"); //getInputValue(wbRange, 3, 2) //Cells(4, 3) // = bs //weibull parameters for stem
      cs = getValueFromNameDbl("i_cs"); //getInputValue(wbRange, 4, 2) //Cells(5, 3) //= cs
                                        //bl = bs: cl = cs //stem and leaf curves equal
      bl = getValueFromNameDbl("i_bl"); //getInputValue(wbRange, 3, 2) //Cells(4, 4) // = bl //weibull parameters for leaf
      cl = getValueFromNameDbl("i_cl"); //getInputValue(wbRange, 4, 2) //Cells(5, 4) // = cl

      lat = getValueFromNameDbl("i_latitude"); //getInputValue(standRange, 1, 1) //Cells(5, 18) //latitude in degree fraction north
      lat = pi / 180 * lat; //converted degrees to radians
      longitude = getValueFromNameDbl("i_longitude"); //getInputValue(standRange, 2, 1) //Cells(6, 18) //longitude in degree fraction west
      tau = getValueFromNameDbl("i_atmTrans"); //Cells(4, 18) //atmospheric transmittance from weather data
      tsncorr = getValueFromNameDbl("i_solarNoon"); //Cells(14, 18) //solar noon correction from weather data
      slope = getValueFromNameDbl("i_slopeI"); //Cells(7, 18) //slope inclination, degrees from horizontal
      slopeasp = getValueFromNameDbl("i_slopeA"); //Cells(8, 18) //slope aspect, counterclockwise degrees from south
      lai = getValueFromNameDbl("i_leafAngleIndex"); //Cells(10, 18) //canopy lai
      xang = getValueFromNameDbl("i_leafAngleParam"); //Cells(9, 18) //leaf angle parameter, CN 15.4

                                                      //'ksatp = Cells(1, 9) //whole plant + rhizosphere KMAX(weibull Kmax) in kg hr - 1 MPa - 1 m - 2 basal area...defined as e / (md - pd)
      gmax = 1000000; //Cells(2, 13) //maximum G, wet soil, vpd = 0, in kg m - 2 hr - 1 basal area
      laperba = getValueFromNameDbl("i_leafPerBasal"); //Cells(1, 13) //initial leaf area per basal area, m2 m - 2
      gmaxl = gmax * (1.0 / laperba) * (1 / 3600.0) * 55.56 * 1000.0; //convert to gmax per leaf area in mmol m-2s-1
                                                                      //'Cells(3, 13) = gmaxl
      alt = getValueFromNameDbl("i_elevation"); //Cells(13, 18) //elevation in m
      patm = 101.325 * pow((1 - 0.0065 * alt / (288.15 + 0.0065 * alt)), 5.257); //atmospheric pressure, T = 15 C, average sealevel patm; approximation
                                                                                 //Call setValueFromName("o_atmP", patm) //Cells(1, 18) = patm //atmospheric pressure in kPa

      leafpercent = getValueFromNameDbl("i_leafPercRes"); //Cells(6, 4) //saturated % of tree R in leaves
                                                          //[HNT] calculate and output the stem and root percent resistances @ ksat -- originally was done by formula in sheet
                                                          // added stemPercent and rootPercent globals at start
      double stempercent = 1.0 / 3.0 * (100.0 - leafpercent);
      double rootpercent = 2.0 / 3.0 * (100.0 - leafpercent);
      //Call setValueFromName("o_stemPercRes", stempercent)
      //Call setValueFromName("o_rootPercRes", rootpercent)

      //ksatp = (leafpercent / 100) * ksatl //kmax of tree
      //Cells(1, 9) = ksatp

      ksatp = getValueFromNameDbl("i_kmaxTree"); //Cells(1, 9) //kmax of tree in kg hr - 1m - 2MPa - 1 per basal area
                                                 //kplantold = ksatp //set kplant initially

                                                 //leafpercent = Cells(6, 4) //% of tree R in leaves
      ksatl = ksatp * (100.0 / leafpercent); //leaf conductance per basal area
      lsc = ksatl * 1.0 / laperba; //lsc in kg hr-1m-2MPa-1//lsc per leaf area in kg hr - 1
                                   //Call setValueFromName("o_leafLSC", lsc / (3600 * 0.00001805)) //Cells(2, 9) = lsc / (3600 * 0.00001805) //lsc converted to mmol s - 1m - 2MPa - 1
      rsatp = 1.0 / ksatp; //convert to resistance
      ksatroot = 1.0 / ((rootpercent / 100.0) * rsatp); //kmax of root system; assumes zero % rhizosphere resistance in WET soil
      ksats = 1.0 / ((stempercent / 100.0) * rsatp); //kmax of stem system
                                                     //Call setValueFromName("o_ksatRoot", ksatroot) //Cells(2, 2) = ksatroot //whole root system kmax
                                                     //Call setValueFromName("o_ksatStem", ksats) //Cells(2, 3) = ksats //stem network kmax
                                                     //Call setValueFromName("o_ksatLeaf", ksatl) //Cells(2, 4) = ksatl //parallel kmax for leaves
                                                     //dbh = ((ksatp * 0.785) / Cells(6, 9)) ^ (1 / (Cells(7, 9) - 2)) //basal diameter in m from k by D allometry
                                                     //height = (99 / Cells(5, 9)) * dbh ^ (2 / 3) //height from H by D allometry and safety factor
                                                     //Cells(3, 9) = height //height in m
                                                     //Cells(3, 7) = dbh
      height = getValueFromNameDbl("i_height"); //Cells(3, 9) //average tree height in m
                                                //dbh = Cells(5, 9) //average tree dbh in m
      pgrav = height * 0.01; //pressure drop from gravity in MPa
      einc = ksatp / 500.0; //e increment in kg hr-1 m-2 basal area for composite curve
                            // aspen pinc = 0.0004;
                            // ponderosa pinc = 0.0007;
                            // original pinc = 0.001;
                            //pinc = 0.0004; // 0.0005; // 0.001; //Cells(8, 13) // MPa increment for global K(P) curves...has to be small enough for NR convergence.
      pinc = getValueFromNameDbl("i_pinc");
      if (pinc <= 0.0)
         pinc = 0.00075; // older parameter sheets don't have pinc, so override if it's zero

      kmin = ksatp / 2000.0; //"instantaneous K" cutoff for global K(P) curves for each element
      aspect = getValueFromNameDbl("i_aspect"); //Cells(5, 12) //max radius of root system per max depth
      layers = getValueFromNameLng("i_layers"); //Cells(7, 12)
                                                //rootfunc = Cells(5, 12) //"y" if function is to be used, "n" if user sets layers
                                                //if rootfunc = "y" { //use the root function
      beta = getValueFromNameDbl("i_rootBeta"); //Cells(6, 12) //root beta for Y = 1 - B ^ d

                                                //for this soil data we want to use the original anchor-offset system

      double layerDepths[20];
      for (k = 1; k <= layers; k++) //set layer depths and % root ksat
      {
         //Cells(8 + k, 2) = 100 * (0.995 / layers) //equal % roots per layer
         //Cells(8 + k, 11) = 0.01 * Log(1 - k * 0.995 / layers) / Log(beta) //lower depth of each layer converted to m
         layerDepths[k] = 0.01 * log(1.0 - k * 0.995 / layers) / log(beta); //lower depth of each layer converted to m
         paramCells[rowLR + k][colLR + 11] = std::to_string(layerDepths[k]);
      }
      //depthmax = Cells(8 + layers, 11) //max depth in meters
      //depthmax = lSheet.Cells(rowLR + layers, colLR + 11) //max depth in meters
      depthmax = layerDepths[layers];
      //calculate transport distances
      //first get vertical distance to biomass center of each layer
      for (k = 1; k <= layers * 2.0; k++)
      {
         vertdistance[k] = 0.01 * log(1 - k * 0.995 / (layers * 2.0)) / log(beta); //get half depths
      }
      i = 0;
      for (k = 1; k <= layers * 2; k += 2) // To layers * 2 Step 2
      {
         i = i + 1;
         vertdistance[i] = vertdistance[k]; //take every other vertdistance
      }
      //now get radial distances
      for (k = 1; k <= layers; k++) //To layers //get thicknesses
      {
         if (k == 1)
         {
            depth[k] = layerDepths[k]; //Cells(8 + k, 11)
         }
         else {
            depth[k] = layerDepths[k] - layerDepths[k - 1]; //lSheet.Cells(rowLR + k, colLR + 11) - lSheet.Cells(rowLR + k - 1, colLR + 11)  //depth is layer thickness in meters
         } // endif
      }
      vol = depth[1] * pi * pow((depthmax * aspect), 2.0); //volume of first, and hence all, layers
                                                           //get radial widths of each layer and transport length
      for (k = 1; k <= layers; k++) //To layers
      {
         radius[k] = pow((vol / (depth[k] * pi)), 0.5); //width in m
         length[k] = radius[k] + vertdistance[k]; //transport distance
         if (k == 1)
            shallow = length[k];
         length[k] = length[k] / shallow; //normalize to shallowest layer
      }
      unknowns = layers + 1; //number of unknowns to be solved for; also dimensions of matrix
      depth[0] = 0; //depth of surface
      rockfrac = getValueFromNameDbl("i_rockFrac"); //Cells(4, 27) //fraction of soil volume as rocks
      rockfrac = 1 - rockfrac; //fraction of volume with no rocks
      for (k = 1; k <= layers; k++) //To layers //read in soil
      {
         //kkmax(k) = Cells(8 + k, 3) //saturated conductivity of soil in kg hr - 1MPa - 1 m - 1
         //a[k] = atof(paramCells[rowLR + k][colLR + 3].c_str());

         a[k] = atof(paramCells[rowLR + k][colLR + 3].c_str()); //van genuchten alpha
         n[k] = atof(paramCells[rowLR + k][colLR + 4].c_str()); //lSheet.Cells(rowLR + k, colLR + 4) //vg n
         kkmax[k] = atof(paramCells[rowLR + k][colLR + 5].c_str());//saturated conductivity of soil in kg hr-1MPa-1 m-1
         thetasat[k] = atof(paramCells[rowLR + k][colLR + 6].c_str()); //theta sat in volume/volume
         thetasat[k] = thetasat[k] * rockfrac; //reduce for actual rock-free fraction of soil
      }
      //now add toplayer (layer 0) of rootless soil 2 cm thick w. same properties as layer 1
      a[0] = a[1];
      n[0] = n[1];
      kkmax[0] = kkmax[1];
      thetasat[0] = thetasat[1];
      depth[0] = 0.02; //sets top layer to 2 cm
                       //now solve for kmax rhizosphere that gives the desired ave % rhizosphere resistance
      rhizotarg = getValueFromNameDbl("i_rhizoPer") / 100.0; //Cells(4, 9) / 100 //average fraction of whole plant resistance in rhizosphere(maximum soil limitation)
      z = 1; //use layer 1 as stand in for whole root system
      ksatr[1] = ksatroot; //set to whole root system
      x = 0.5; //start by finding kmaxrh at 0.5 MPa...a deliberate under-shoot
      rootr = 1.0 / wbr(x);
      rstem = 1.0 / wbs(x);
      rleaf = 1.0 / wbl(x);
      rplant = rootr + rstem + rleaf; //rplant here is just the xylem part
      rhizor = rplant * (rhizotarg / (1.0 - rhizotarg)); //solve for what rhizor has to be at the target
      vp = 1.0 / (pow((a[z] * x), n[z]) + 1); //van genuchten terms // vp = 1 / ((a(z) * x) ^ n(z) + 1) 
      vgterm = pow(vp, ((n[z] - 1) / (2.0 * n[z]))) * pow((pow((1 - vp), ((n[z] - 1) / n[z])) - 1), 2.0); //van genuchten terms // vgterm = vp ^ ((n[z] - 1) / (2 * n[z])) * ((1 - vp) ^ ((n[z] - 1) / n[z]) - 1) ^ 2
      kmaxrh[1] = (1.0 / rhizor) / vgterm; //solve for kmaxrh[1]
      kinc = kmaxrh[1] * 0.1;

      //std::cout << "do 1 begin, kinc = " << kinc << " rhizotarg = " << rhizotarg * 100.0 << " pinc = " << pinc << std::endl;

      do //loop through rhizosphere kmax
      {
         kmaxrh[1] = kmaxrh[1] + kinc; //increase from deliberate undershoot
         x = 0; //
         sum = 0;
         do //loop through pressures
         {
            x = x + 0.1;
            rootr = 1.0 / wbr(x);
            rstem = 1.0 / wbs(x);
            rleaf = 1.0 / wbl(x);
            rhizor = 1.0 / vg(x);
            rplant = rootr + rstem + rleaf + rhizor;
            rrhizofrac = rhizor / rplant; //fraction of resistance in rhizosphere
            sum = sum + rrhizofrac; //add up fractions
         } while (!(1.0 / rplant < kmin)); //Loop Until 1 / rplant < kmin //average over full range
         sum = sum / (x / 0.1); //average fraction
      } while (!(sum < rhizotarg)); // Until sum < rhizotarg //loop until desired soil limitation is reached

                                    //std::cout << "do 1." << std::endl;

      kmaxrh[1] = kmaxrh[1] / layers; //divide whole root rhizokmax into equal portions for each layer
                                      //end of soil limitation adjustment
                                      //now set soil layer parameters based on aroot of entire root system
      for (z = 1; z <= layers; z++) //z = 1 To layers
      {
         kmaxrh[z] = kmaxrh[1]; //soil to root MAXIMUM conductance in kg hr-1 MPa-1//re - set for individual layers
                                //lSheet.Cells(rowLR + z, colLR + 7) = kmaxrh(z) //note: this is kMAX...at P=0; not at saturated PD
         paramCells[rowLR + z][colLR + 7] = std::to_string(kmaxrh[z]); //note: this is kMAX...at P=0; not at saturated PD
      }
      t = 0;
      //loop to find ksatroot for each layer
      coef = 0.0;
      do
      {
         coef = coef + 0.01;
         sum = 0.0;
         for (k = 1; k <= layers; k++)//k = 1 To layers //soil layers from top to bottom
         {
            ksatr[k] = coef / length[k]; //assumes ksatr proportional to biomass/length
            sum = sum + ksatr[k];
         }
      } while (!(sum > ksatroot));
      //Loop Until sum > ksatroot //loop until each layer adds to total

      //std::cout << "do 2." << std::endl;

      for (k = 1; k <= layers; k++)//k = 1 To layers
      {
         //lSheet.Cells(rowLR + k, colLR + 2) = 100 * (ksatr(k) / ksatroot);
         //lSheet.Cells(rowLR + k, colLR + 9) = ksatr(k);
         paramCells[rowLR + k][colLR + 2] = std::to_string(100.0 * (ksatr[k] / ksatroot));
         paramCells[rowLR + k][colLR + 9] = std::to_string(ksatr[k]);
      }
      refilling = getValueFromNameStr("i_refilling"); //Cells(4, 12) //"y/n"
                                                      //A-ci curve parameters...A = aaa(1-bbb^ci)
                                                      //par = Cells(3, 21)//now read in each time step
      qmax = getValueFromNameDbl("i_qMax"); //Cells(3, 21) //quantum yield of electron transport, moles e per mols photons
      ca = getValueFromNameDbl("i_co2AmbPPM"); //Cells(2, 21) //ambient co2 in ppm
      ca = ca * 0.000001; //ambient co2 in moles per mole
                          //Cells(1, 21) = ca * patm * 1000 //ambient co2 in Pa
                          //Call setValueFromName("o_co2AmbPa", ca * patm * 1000) // this is a pure output from the Excel version, OK to disable here
                          // the value in this cell will not be referenced again by the model
      vmax25 = getValueFromNameDbl("i_vmax25"); //Cells(4, 21) //vmax at 25C
                                                //vmax25initial = vmax25
      jmax25 = getValueFromNameDbl("i_jmax25"); //Cells(5, 21) //jmax at 25C
                                                //rday25 = Cells(8, 21) //day respiration at 25C, = 0.015vmax cf.collatz et al.
      kc25 = getValueFromNameDbl("i_kc25"); //Cells(6, 21) //m - m constant for CO2 in mole fraction at 25C
      ko25 = getValueFromNameDbl("i_ko25"); //Cells(7, 21) //m - m constant for O2 in mole fraction at 25C
      comp25 = getValueFromNameDbl("i_comp25"); //Cells(8, 21) //photorespiratory compensation point in mole fraction at 25C
      thetac = getValueFromNameDbl("i_thetaC"); //Cells(9, 21) //shape factor for A - ci colimitation, 0.98
                                                //airtemp = Cells(3, 18) //in C...NOW INPUTED PER TIME STEP
                                                //maxvpd = (101.3 / patm) * (-0.0043 + 0.01 * Exp(0.0511 * airtemp)) //saturated mole fraction of vapor in air
                                                //Cells(15, 1) = maxvpd * patm //maximum possible vpd
                                                //rabs = Cells(4, 18) //long and short wave radiation absorbed PER TIME STEP
                                                //wind = Cells(5, 18) //wind speed PER TIME STEP
      leafwidth = getValueFromNameDbl("i_leafWidth") * 0.72; //Cells(2, 18) * 0.72 //leaf width x factor = characteristic dimension(campbell and norman)
      emiss = getValueFromNameDbl("i_emiss"); //Cells(3, 18) //long wave emissivity
      havmax = getValueFromNameDbl("i_havmax"); //Cells(10, 21) //these are all temp - dependency parameters from Leunig 2002
      hdvmax = getValueFromNameDbl("i_hdvmax"); //Cells(11, 21)
      svvmax = getValueFromNameDbl("i_svvmax"); //Cells(12, 21)
      hajmax = getValueFromNameDbl("i_hajmax"); //Cells(10, 23)
      hdjmax = getValueFromNameDbl("i_hdjmax"); //Cells(11, 23)
      svjmax = getValueFromNameDbl("i_svjmax"); //Cells(12, 23)
      lightcurv = getValueFromNameDbl("i_lightCurv"); //Cells(13, 21)
      baperga = getValueFromNameDbl("i_baperga") * 0.0001; //Cells(3, 27) * 0.0001 //basal area per ground area converted from m2 / Ha to m2 / m2
      ffc = getValueFromNameDbl("i_fieldCapPercInit") / 100.0; //Cells(2, 27) / 100 //fraction of field capacity for starting the season
      fieldcapfrac = getValueFromNameDbl("i_fieldCapFrac"); //Cells(1, 27) //fraction that field capacity is of saturation(minus residual)
                                                            //timestep = Cells(4, 27) //timestep in fraction of hour
      pground = getValueFromNameDbl("i_gWaterP"); //Cells(5, 27) //ground water pressure
      grounddistance = getValueFromNameDbl("i_gWaterDist"); //Cells(6, 27) //distance to ground water source in m
      ground = getValueFromNameStr("i_gWaterEnable"); //Cells(7, 27) //groundwater flow, yes or no
      soilred = getValueFromNameStr("i_soilRedEnable"); //Cells(8, 27) //turns off / on soil redistribution routine
      sevap = getValueFromNameStr("i_soilEvapEnable"); //Cells(9, 27) //turns off / on soil evaporation routine
      soilabssol = getValueFromNameDbl("i_soilAbsSol"); //Cells(11, 18) //absorptivity of soil surface for solar
      rough = 0.01; //soil Zm, eqn 14.9, using table 5.1 for smooth surface, cm
      zdispl = 6.5 * rough; //soil d, eqn 14.9, using d = 6.5 Zm, eq 5.2,5.3
      xh = getValueFromNameDbl("i_soilXHeight") * 100.0; //Cells(12, 18) * 100 //height above soil surface for understory wind and gh in cm
      zh = 0.2 * rough; //roughness for temperature
                        //rainsim = getValueFromName("i_rainEnable"); //Cells(11, 27) //y/n...turns on simulated rain
      raining = getValueFromNameStr("i_rainEnable"); //[HNT] enable or disable OBSERVED rainfall (if //n// will ignore rain data)
      if (raining == "n") {
         rainEnabled = false;
      }
      else {
         rainEnabled = true; //enabled is made the "else" case becasue the "default" state is to process rain, so any input OTHER than //n// enables
      } // endif

        // opposite arrangement here: set GS data usage FALSE unless it's explicitly == "y"
        // because this requires an extra data file we may not have
      if (stage_ID == STAGE_ID_NONE || stage_ID == STAGE_ID_HIST_STRESS || stage_ID == STAGE_ID_FUT_STRESS || stage_ID == STAGE_ID_FUT_STRESS_NOACCLIM)
      {
         if (getValueFromNameStr("i_useGSDataStress") == "y") {
            useGSData = true;
         }
         else {
            useGSData = false;
         } // endif
      }
      else if (stage_ID == STAGE_ID_HIST_OPT || stage_ID == STAGE_ID_FUT_OPT)
      {
         if (getValueFromNameStr("i_useGSDataOpt") == "y") { // different parameter name
            useGSData = true;
         }
         else {
            useGSData = false;
         } // endif
      }
      else // impossible unknown stage failsafe
         useGSData = false;
   }

   short initModelVars()
   {
      // this can be called on the start of every iteration BUT NOT ON NEW YEARS
      isNewYear = true;
      // this is a good place for general initialization too
      gs_yearIndex = 0;
      gs_prevDay = 0;
      gs_inGrowSeason = false;
      gs_doneFirstDay = false;
      year_cur = 0;
      year_start = 0;
      yearVal = 0;

      return -1;
   }

   short cleanModelVars()
   {
      // model program vars
      long iii = 0;
      long qqq = 0;
      e = 0;
      br = 0;
      cr = 0;

      x = 0;

      ksh = 0;

      for (iii = 0; iii < 6; iii++)
      {
         //6 arrays 1-d
         n[iii] = 0;
         a[iii] = 0;
         ksatrh[iii] = 0;
         kmaxrh[iii] = 0;
         ksatr[iii] = 0;
         kmaxr[iii] = 0;
         depth[iii] = 0;
         pcritr[iii] = 0;
         prh[iii] = 0;
         kminroot[iii] = 0;
         length[iii] = 0;
         radius[iii] = 0;
         layer[iii] = 0;
         tlayer[iii] = 0;
         layerfailure[iii] = "ok";
         tlayerfailure[iii] = "ok";

         soilf[iii] = 0;
         kkmax[iii] = 0;
         thetasat[iii] = 0;
         swclimit[iii] = 0;
         if (!(useGSData && gs_yearIndex > 0))
         {
            thetafc[iii] = 0;
            thetafracres[iii] = 0;
            thetafracfc[iii] = 0;

            water[iii] = 0;
            fc[iii] = 0;
         }

         pcritrh[iii] = 0;
         pd[iii] = 0;
         soilredist[iii] = 0;

         for (qqq = 0; qqq < 1000001; qqq++)
         {
            // 1-D 100k and 1mil arrays
            if (iii == 0)
            {
               // 100k arrays
               if (qqq < 100001)
               {
                  psynshmd[qqq] = 0;
                  psynmd[qqq] = 0;
                  psynsh[qqq] = 0;
                  amaxfracsh[qqq] = 0;
                  leaftempsh[qqq] = 0;
                  lavpdsh[qqq] = 0;
                  rdaysh[qqq] = 0;
                  gcanwsh[qqq] = 0;
                  gcancsh[qqq] = 0;
                  cinsh[qqq] = 0;
                  rday[qqq] = 0;
                  leaftemp[qqq] = 0;
                  lavpd[qqq] = 0;
                  amaxfrac[qqq] = 0;
                  dpa[qqq] = 0;
                  psyn[qqq] = 0;
                  gcanw[qqq] = 0;
                  gcanc[qqq] = 0;
                  dedp[qqq] = 0;
                  pleaf[qqq] = 0;
                  kplant[qqq] = 0;
                  eplant[qqq] = 0;
                  kstem[qqq] = 0;
                  pstem[qqq] = 0;
                  proot[qqq] = 0;
                  kleaf[qqq] = 0;
                  tel[qqq] = 0;
                  tes[qqq] = 0;
                  el[qqq] = 0;
                  es[qqq] = 0;
               }
               // 1 mil arrays
               kloss[qqq] = 0;
               klossv[qqq] = 0;
               pleafv[qqq] = 0;
               eplantl[qqq] = 0;
               dedpf[qqq] = 0;
               cin[qqq] = 0;
            }
            // 6,100k 2d arrays 
            if (qqq < 100001) //don't mix up the 100k and 1mil arrays
            {
               er[iii][qqq] = 0;
               kr[iii][qqq] = 0;
               erh[iii][qqq] = 0;

               krh[iii][qqq] = 0;
               tkr[iii][qqq] = 0;
               ter[iii][qqq] = 0;
               elayer[iii][qqq] = 0;
               prhizo[iii][qqq] = 0;
               kroot[iii][qqq] = 0;
            }
            // 6,1 mil arrays 2d
         }
      }

      // 7 arrays
      for (iii = 0; iii < 7; iii++)
      {
         indx[iii] = 0;
         func[iii] = 0;
         vv[iii] = 0;
         dfrhdprh[iii] = 0;
         dfrdprh[iii] = 0;
         dfrhdpr[iii] = 0;
      }

      // lazy memsets for anything that doesn't fit -- usually good enough
      std::memset(jmatrix, 0, sizeof(jmatrix));
      memset(vertdistance, 0, sizeof(vertdistance));
      memset(inp, 0, sizeof(inp));
      //jmatrix[7][7] = 0;
      //vertdistance[11] = 0;
      //inp[101];

      bs = 0;
      cs = 0;
      bl = 0;
      cl = 0;
      ksats = 0;
      ksatl = 0;

      pcrits = 0;
      pcritl = 0;
      ps = 0;

      pinc = 0;
      einc = 0;
      ksatp = 0;
      rsatp = 0;
      ksatroot = 0;
      vp = 0;
      s = 0;
      del = 0;
      eps = 0;
      epsx = 0;
      sthresh = 0;
      sum = 0;
      sums = 0;
      olds = 0;
      kmin = 0;

      elow = 0;
      ehigh = 0;
      plow = 0;

      estart = 0;
      efinish = 0;
      flow = 0;
      dfrdpr = 0;
      frt = 0;
      klow = 0;
      khigh = 0;
      kupper = 0;
      klower = 0;
      pr = 0;
      threshold = 0;
      kfactor;
      aamax = 0;
      dum = 0;
      kminleaf = 0;
      kminstem = 0;
      md = 0;
      ecritsystem = 0;
      pcritsystem;
      phigh = 0;
      failspot = "";
      setting = "";
      refilling = "";
      ground = "";
      soilred = "";
      dp = 0;
      kminplant = 0;
      toplayer;
      p1 = 0;
      p2 = 0;
      pl = 0;
      plold = 0;
      predawn = 0;
      dedpl = 0;
      prtarget = 0;
      kmaxs = 0;
      pstarget = 0;
      kmaxl = 0;
      pltarget = 0;
      rhizsat = 0;
      rootsat = 0;
      sum2 = 0;
      beta = 0;
      rhizok = 0;
      patm = 0;
      vpd = 0;
      dedplzero;
      gcmd = 0;
      dpmax = 0;
      gmax = 0;
      vpdsat = 0;
      pstop = 0;
      halt = 0;
      failure = 0;
      k = 0;
      check = 0;
      t = 0;
      j = 0;
      it = 0;
      f = 0;
      //dd = 0;
      tmax = 0;
      tnm = 0;
      test = 0;
      layers = 0;
      d = 0;
      z = 0;
      p = 0;
      pmax = 0;
      unknowns = 0;
      i = 0;
      imax = 0;
      ii = 0;
      ll = 0;
      total = 0;
      kplantold = 0;
      frac = 0;
      eincdef = 0;
      aspect = 0;
      vol = 0;
      depthmax = 0;
      shallow = 0;
      coef = 0;
      rhizotarg = 0;
      rplant = 0;
      rstem = 0;
      rootr = 0;
      rleaf = 0;
      rhizor = 0;
      rrhizofrac = 0;
      kinc = 0;
      kplantmax = 0;
      vgterm = 0;
      ci = 0;
      ca = 0;
      marker = 0;
      var = 0;
      transpiration = 0;
      dedpmax = 0;
      lavpdmd = 0;
      psynact = 0;
      psynmax = 0;
      dpamax;
      grad = 0;
      gha = 0;
      numerator = 0;
      denominator = 0;
      lambda = 0;
      emiss = 0;
      airtemp = 0;
      rabs = 0;
      laperba = 0;
      par = 0;
      qmax = 0;
      vmax25 = 0;
      kc25 = 0;
      ko25 = 0;
      comp25 = 0;
      theta = 0;
      wind = 0;
      leafwidth = 0;
      comp = 0;
      vmax = 0;
      kc = 0;
      ko = 0;
      je = 0;
      jc = 0;
      lavpdc1 = 0;
      lavpdsum = 0;
      lavpdh = 0;
      gmaxl = 0;
      cinc = 0;
      jmax = 0;
      jmax25 = 0;
      havmax = 0;
      hdvmax = 0;
      svvmax = 0;
      hajmax = 0;
      hdjmax = 0;
      svjmax = 0;
      jact = 0;
      maxvpd = 0;
      maxkloss = 0;
      chalk = 0;
      skip = 0;
      psynstop = 0;
      totalv = 0;
      lsc = 0;
      timestep = 0;
      leafpercent = 0;
      dbh = 0;
      height = 0;
      pgrav = 0;
      ffc = 0;
      pground = 0;
      grounddistance = 0;
      baperga = 0;
      layerflow = 0;
      groundwater = 0;
      deficit = 0;
      rain = 0;
      sumrain = 0;
      store = 0;
      pend = 0;
      groundflow = 0;
      rday25 = 0;

      if (!(useGSData && gs_yearIndex > 0))
      {
         waterold = 0;
         waternew = 0;
         waterchange = 0;
      }

      drainage = 0;
      gwflow = 0;
      thetac = 0;
      lat = 0;
      longitude = 0;
      slope = 0;
      slopeasp = 0;
      lai = 0;
      xang = 0;
      fet = 0;
      et = 0;
      sm = 0;
      lc = 0;
      tsn = 0;
      sindec = 0;
      dec = 0;
      cosdec = 0;
      tim = 0;
      tod = 0;
      coszen = 0;
      zen = 0;
      cosaz = 0;
      az = 0;
      m = 0;
      sp = 0;
      sb = 0;
      sd = 0;
      st = 0;
      cloud = 0;
      obssolar = 0;
      fcd = 0;
      kbe = 0;
      kbezero = 0;
      mleafang = 0;
      rad = 0;
      k1 = 0;
      t1 = 0;
      told = 0;
      t2 = 0;
      kd = 0;
      qd = 0;
      qds = 0;
      qdt = 0;
      qb = 0;
      qbt = 0;
      qsc = 0;
      qsh = 0;
      qsl = 0;
      laisl = 0;
      laish = 0;
      parsh = 0;
      parsl = 0;
      parbottom = 0;
      nirsh = 0;
      nirsl = 0;
      sshade = 0;
      ssun = 0;
      sbottom = 0;
      ssunb = 0;
      ssund = 0;
      sref = 0;
      ppfd = 0;
      ea = 0;
      eac = 0;
      la = 0;
      lg = 0;
      jd = 0;
      idiot = 0;
      o = 0;
      mdsh = 0;
      cincsh = 0;
      lavpdmdsh = 0;
      gcmdsh = 0;
      psynactsh = 0;
      transpirationsh = 0;
      haltsh = 0;
      psynmaxsh = 0;
      transpirationtree = 0;
      anetsh = 0;
      anettree = 0;
      atree = 0;
      anet = 0;
      rabssoil = 0;
      alt = 0;
      tsncorr = 0;
      rha = 0;
      rhs = 0;
      soilevap = 0;
      soilabssol = 0;
      rough = 0;
      zdispl = 0;
      xh = 0;
      zh = 0;
      us = 0;
      mdensair = 0;
      soilep = 0;
      emission = 0;
      soiltemp = 0;
      night = "";
      sevap = "";
      pet = "";
      rainsim = "";
      nightBool = 0;
      tau = 0;
      minwind = 0;
      runmean = 0;
      weird = 0;
      sign = 0;
      xx = 0;
      runs = 0;
      ticks = 0;
      rmean = 0;
      cutoff = 0;
      dedplmin = 0;
      amaxmax = 0;
      lightcurv = 0;
      emd = 0;
      leaftmd = 0;
      leaftshmd = 0;
      lavpdshmd = 0;
      gcanwmd = 0;
      gcancmd = 0;
      rdaymd = 0;
      psynmaxmd = 0;
      cinmd = 0;
      gcanwshmd = 0;
      gcancshmd = 0;
      rdayshmd = 0;
      cinshmd = 0;
      psynmaxshmd = 0;
      prinitial = 0;
      dpasun = 0;
      lightcomp = 0;
      fieldcapfrac = 0;
      rockfrac = 0;
      initialthreshold = 0;
      kref = 0;
      lscref = 0;
      pdref = 0;
      leafpercentref = 0;
      kpday1 = 0;
      kxday1 = 0;
      runoff = 0;
      kmaxits = 0;
      bst = 0;
      startbst = 0;
      scenario = 0;
      kmaxset = 0;
      leafpercentset = 0;
      return -1; // true in VBA langauge
   }

   short setIterationParams(long &iter)
   {
      if (iter < 1000 && iter > -1)
      {
         //gs_yearIndex = iter;
         iter_Counter = iter;
         return -1; // true in short language
      }
      else
      {
         return 0;
      }
   }

   bool isInGrowSeasonSimple()
   {
      if (useGSData)
      {
         if (gs_yearIndex >= 0 && gs_yearIndex < 100)
         {
            if (gs_ar_years[gs_yearIndex] > 0)
            {
               if (jd >= gs_ar_starts[gs_yearIndex] && jd <= gs_ar_ends[gs_yearIndex])
               {
                  return true;
               }
            }
         }
         return false;
      }
      else
      {
         return true; // if the GS limits are disabled, we're always in the growing season
      }
   }

   void updatecurves() //'resets element E(P) curves
   {
      //'if k<kmin, re-assign e//'s on element curve by back-calculating
      for (z = 1; z <= layers; z++)//z = 1 To layers
      {
         if (true)
         {
            if (kroot[z][halt] < kminroot[z])
            {
               kminroot[z] = kroot[z][halt];
               phigh = int(proot[halt] / pinc) + 1; //'pressure datum just above the target
               for (k = phigh; k >= 0; k--)//k = phigh To 0 Step -1 //'back-calculate e//'s
               {
                  er[z][k] = er[z][phigh] - kminroot[z] * pinc * (phigh - k);
                  kr[z][k] = kminroot[z]; //'back-calculate KR(Z,K) too for roots (not stem or leaves)
               } //EndFor  k
            } //EndIf//
         }
      } //EndFor  z
      if (kstem[halt] < kminstem)
      {
         kminstem = kstem[halt];
         phigh = int(pstem[halt] / pinc) + 1;
         for (k = phigh; k >= 0; k--)//k = phigh To 0 Step -1 //'back-calculate e//'s
         {
            es[k] = es[phigh] - kminstem * pinc * (phigh - k);
         } //EndFor  k
      } //EndIf//
      if (kleaf[halt] < kminleaf)
      {
         kminleaf = kleaf[halt];
         phigh = int(pleaf[halt] / pinc) + 1;
         for (k = phigh; k >= 0; k--)//k = phigh To 0 Step -1 //'back-calculate e//'s
         {
            el[k] = el[phigh] - kminleaf * pinc * (phigh - k);
         } //EndFor  k
      } //EndIf//
        //'if kplant[halt] < kminplant Then kminplant = kplant[halt]NOTE: kplant CAN go up because of rhizosphere recovery!
   }

   void solarcalc() //'gets radiative terms for energy balance and assimilation
   {
      //'j = Cells(14 + i, 3) //'julian day
      fet = 279.575 + 0.9856 * jd; //'jd is julian day, fet is factor for CN eqn 11.4
      fet = fet * pi / 180.0; //'convert to radians
      et = (-104.7 * sin(fet) + 596.2 * sin(2 * fet) + 4.3 * sin(3 * fet) - 12.7 * sin(4 * fet) - 429.3 * cos(fet) - 2 * cos(2 * fet) + 19.3 * cos(3 * fet)) / 3600.0; //'"equation of time" in fraction of HOURS C&N 11.4
                                                                                                                                                                       //'long = Cells(11, 4) //'longitude in degree fraction W
      sm = 15 * int(longitude / 15.0); //'standard meridian east of longitude
                                       //'lc = 0.0666667 * (sm - longitude) //'CN 11.3 for longitude correction in fraction of hours
      tsn = 12 - tsncorr - et; //'time of solar noon in hour fraction from midnight, CN 11.3
                               //'Cells(14 + i, 5) = tsn //'output solar noon
      sindec = pi / 180.0 * (356.6 + 0.9856 * jd);
      sindec = sin(sindec);
      sindec = pi / 180.0 * (278.97 + 0.9856 * jd + 1.9165 * sindec);
      sindec = 0.39785 * sin(sindec); //'sine of the angle of solar declination
                                      //'lat = Cells(10, 4) //'latitude, fraction of degrees N
      dec = atan(sindec / pow((-sindec * sindec + 1), 0.5)); //'arcsin of sindec in radians, dec is solar declination
                                                             //'Cells(14 + i, 6) = dec * 180 / pi //'output solar declination
      cosdec = cos(dec);
      //'timst = Cells(14 + i, 4) //'local standard time, hour fraction from midnight
      tim = 15 * (tod - tsn);
      tim = pi / 180.0 * tim; //'convert to radians
      coszen = sin(lat) * sindec + cos(lat) * cosdec * cos(tim); //'cos of zenith angle of sun from overhead, CN11.1
      zen = atan(-coszen / pow((-coszen * coszen + 1), 0.5)) + 2 * atan(1); //'zenith in radians
                                                                            //'if zen < 1.57 { Cells(14 + i, 7) = zen * 180 / pi //'output when sun//'s up (zen<90)
      cosaz = -(sindec - cos(zen) * sin(lat)) / (cos(lat) * sin(zen)); //'cos of azimuth angle measured counterclockwize from due south CN 11.5
      if (cosaz < -1)
         cosaz = -1; //'keep it in limits
      if (cosaz > 1)
         cosaz = 1; //'ditto
      if (cosaz == 1 || cosaz == -1) { //'keeps stupid acos eqn from crashing
         if (cosaz == 1)
            az = 0;
         if (cosaz == -1)
            az = 3.14159; //'180 in radians
      }
      else {
         az = atan(-cosaz / pow((-cosaz * cosaz + 1), 0.5)) + 2 * atan(1); //'solar az in radians
      } //Endif//
      if (tod > tsn)
         az = 6.28319 - az; //'correct for afternoon hours!
                            //'if zen < 1.57 { Cells(14 + i, 8) = az * 180 / pi//'output azimuth during day
                            //'dayl = (-Sin(lat) * sindec) / (Cos(lat) * Cos(dec))
                            //'dayl = Atn(-dayl / Sqr(-dayl * dayl + 1)) + 2 * Atn(1)
                            //'dayl = 2 * dayl / 15
                            //'Cells(14 + i, 9) = dayl * 180 / pi
      if (zen * 180 / pi < 90) { //'sun//'s up: calculate solar radiation
                                 //'pa = Cells(12, 4) //'atmospheric p in kpa
                                 //'Cells(16 + dd, 8) = zen * 180 / pi //'zenith angle in degrees
                                 //'Cells(16 + dd, 9) = lat
                                 //'night = "n" //'its officially day
         m = patm / (101.3 * cos(zen)); //'CN 11.12
                                        //'spo = Cells(9, 4) //'solar constant
                                        //'tau = Cells(8, 4) //'transmittance, CN 11.11
         sp = solar * pow(tau, m); //'direct beam irradiance Wm-2
         sb = sp * cos(zen); //'direct beam irradiance on horizontal surface
                             //'Cells(14 + i, 11) = sb
         sd = 0.3 * (1 - pow(tau, m)) * solar * cos(zen); //'clear sky diffuse radiation
                                                          //'Cells(14 + i, 12) = sd
         st = sd + sb; //'total horizontal irradiance from sun (w/o reflected radiation)
         cloud = solar * pow(0.4, m) * cos(zen); //'overcast threshold
                                                 //'stobs = Cells(13, 4) //'observed solar radiation on the horizontal Wm-2
         if (obssolar > 0) { //'we//'ve got solar data
            if (obssolar < st) { //'we//'ve got clouds
               if (obssolar > cloud) { //'we//'ve got partial clouds
                  fcd = 1 - (obssolar - cloud) / (st - cloud); //'fraction for converting beam to diffuse
                  sd = sd / st + fcd * sb / st; //'diffuse/total rato
                  sd = sd * obssolar; //'multiply ratio by total observed to get total diffuse
                  sb = obssolar - sd; //'leftover beam
                  st = obssolar; //'reset to stobs
               }
               else { //'its all clouds
                  sd = obssolar;
                  sb = 0;
                  st = obssolar;
               } //Endif//
            } //Endif// //'if no clouds, everything//'s already set
         } //Endif// //'if no solar data, we assume no clouds
           //'calculate reflected light as if it is equal to light at bottom of canopy
           //'xang = Cells(12, 11) //'leaf angle parameter, CN 15.4
           //'if zen < 1.57 { //'sun//'s up:
         double evenMoreTest = tan(zen);
         evenMoreTest = pow((tan(zen)), 2.0);
         evenMoreTest = pow(xang, 2.0);
         evenMoreTest = pow((pow(xang, 2.0) + pow((tan(zen)), 2.0)), 2.0);
         double testNum = pow((pow(xang, 2.0) + pow((tan(zen)), 2.0)), 2.0);
         double testDenom = (xang + 1.774 * pow((xang + 1.182), -0.733));
         kbe = pow((pow(xang, 2.0) + pow((tan(zen)), 2.0)), 0.5) / (xang + 1.774 * pow((xang + 1.182), -0.733)); //'beam extinction coefficient CN 15.4
                                                                                                                 //kbe = Sqr(xang ^ 2 + (Tan(zen)) ^ 2) / (xang + 1.774 * (xang + 1.182) ^ -0.733) 'beam extinction coefficient CN 15.4
         kbezero = xang / (xang + 1.774 * pow((xang + 1.182), -0.733)); //'beam extinction for zen=0(overhead
         mleafang = atan(-kbezero / pow((-kbezero * kbezero + 1), 0.5)) + 2 * atan(1); //'mean leaf angle in radians
                                                                                       //'Cells(14 + i, 20) = kbe
                                                                                       //'Cells(14 + i, 21) = mleafang * 180 / pi //'mean leaf angle in degrees
                                                                                       //'lai = Cells(13, 11) //'canopy leaf area index
                                                                                       //'abspar = Cells(7, 17) //'absorptivity for PAR of leaves
                                                                                       //'abssol = Cells(8, 17) //'absorptivity for total solar of leaves
                                                                                       //'gdbeamf = Exp(-Sqr(abssol) * kbe * lai) //'CN 15.6, fraction of solar beam radiation reaching ground
                                                                                       //'Cells(14 + i, 22) = gdbeamf
                                                                                       //'now solve for kd (diffuse extinction) by integrating beam over all possible zenith angles from 0 to 90, CN 15.5
         rad = 0;
         sum = 0;
         k1 = pow((pow(xang, 2) + pow((tan(rad)), 2)), 0.5) / (xang + 1.774 * pow((xang + 1.182), -0.733));  //'beam extinction coefficient CN 15.4
         t1 = exp(-k1 * lai); //'transmittance CN 15.1
         told = t1 * sin(rad) * cos(rad); //'integral function
         do
         {
            rad = rad + 0.015708; //'0.9 degree intervals, zenith angle in radians
            k1 = pow((pow(xang, 2) + pow((tan(rad)), 2)), 0.5) / (xang + 1.774 * pow((xang + 1.182), -0.733));  //'beam extinction coefficient CN 15.4
            t1 = exp(-k1 * lai); //'transmittance CN 15.1
            t2 = t1 * sin(rad) * cos(rad); //'integral function
            sum = sum + (t2 + told) / 2.0 * 0.015708; //'integral sum
            told = t2; //'reset
         } while (!(rad > 1.5708));
         //Loop Until rad > 1.5708 //'loop until 90 degrees
         sum = sum * 2; //'complete summing
         kd = -log(sum) / lai; //'extinction coefficient for diffuse radiation
                               //'Cells(14 + i, 23) = kd //'output
                               //'now...compute shaded leaf ppfd, q denotes a ppfd
         qd = 0.45 * sd * 4.6; //'converts total solar diffuse to PPFD diffuse in umol m-2 s-1
         qds = qd * (1 - exp(-(pow(abspar, 0.5) * kd * lai))) / (pow(abspar, 0.5) * kd * lai); //'mean diffuse irradiance for shaded leaves CN p. 261
         qdt = qd * exp(-(pow(abspar, 0.5) * kd * lai)); //'diffuse irradiance at bottom of canopy, CN p. 255, Eqn 15.6
         qb = 0.45 * sb * 4.6; //'converts total solar beam to PPFD
         qbt = qb * exp(-(pow(abspar, 0.5) * kbe * lai)); //'direct AND downscattered PPFD at bottom of canopy, CN 15.6
         qb = qb * exp(-(kbe * lai)); //'direct PPFD at bottom of canopy,CN 15.1
         qsc = (qbt - qb) / 2.0; //'average backscattered beam on shaded leaves
         qsh = qds + qsc; //'average PPFD incident (not absorbed!) on shaded leaves
         qb = 0.45 * sb * 4.6; //'re-set qb to top of canopy
                               //'now get sunlit leaf ppfd
         qsl = kbe * qb + qsh; //'average PPFD incident (not absorbed) on sunlit leaves
                               //'now get sunlit vs. shaded lai
         laisl = (1 - exp(-(kbe * lai))) / kbe; //'sunlit lai
         laish = lai - laisl; //'shaded lai
                              //'now get sun and shade PAR, NIR, & longwave in preparation for energy balance
                              //'in par range
         parsh = qsh / 4.6; //'incident PAR, Wm-2, shaded leaves, 100% diffuse
         parsl = qsl / 4.6;//'incident PAR, sunlit leaves
         parbottom = (qbt + qdt) / 4.6; //'PAR making it through the canopy
                                        //'in near-infrared range (assume same equations as for PAR, but different absorptances and incoming fluxes in Wm-2
                                        //'absnir = Cells(9, 17) //'absorptivity of leaves to NIR
         qd = 0.55 * sd; //'converts total solar diffuse to NIR diffuse in Wm-2
         qds = qd * (1 - exp(-(pow(absnir, 0.5) * kd * lai))) / (pow(absnir, 0.5) * kd * lai); //'mean diffuse nir irradiance for shaded leaves CN p. 261
         qdt = qd * exp(-(pow(absnir, 0.5) * kd * lai)); //'diffuse NIR irradiance at bottom of canopy, CN p. 255, Eqn 15.6
         qb = 0.55 * sb; //'converts total solar beam to NIR
                         //'qb = 1600
         qbt = qb * exp(-(pow(absnir, 0.5) * kbe * lai)); //'direct AND downscattered NIR at bottom of canopy, CN 15.6
         qb = qb * exp(-(kbe * lai)); //'direct NIR at bottom of canopy,CN 15.1, using same extinction coefficient as for PAR
         qsc = (qbt - qb) / 2.0; //'average backscattered beam on shaded leaves
         nirsh = (qds + qsc); //'incident NIR on shaded leaves, 100% diffuse
         qb = 0.55 * sb; //'re-set nir qb to top of canopy to get sunlit leaves
         nirsl = kbe * qb + nirsh; //'average incident NIR on sunlit leaves
         sshade = parsh + nirsh; //'total solar incident on shaded leaves, 100% diffuse
         ssun = parsl + nirsl; //'total solar incident on sunlit leaves
         sbottom = parbottom + qdt + qbt; //'total solar at bottom of canopy
         ssunb = sb / st * ssun; //'beam solar on sunlit (an approximation)
         ssund = sd / st * ssun;//'diffuse solar on sunlit (approximation)
                                //'abssolar = Cells(8, 17) //'absorptivity of leaves for total solar (0.5)
         sref = (1 - absolar) * sshade; //'reflected light...this used for sun/shade dichotomy
         sref = (1 - absolar) * sbottom; //'reflected light for monolayer version
                                         //'these below are used for monolayer version:
         par = 0.45 * st; //'wm-2 in par wavelength...45% of total solar
         ppfd = par * 4.6; //'assumes 4.6 moles photons per Joule conversion factor
      }
      else { //'sun//'s down
         sp = 0; sb = 0; sd = 0; st = 0; sref = 0; par = 0; ppfd = 0; sref = 0; ssun = 0; sshade = 0;
         qsh = 0; qsl = 0; ssunb = 0; ssund = 0; laisl = 0; laish = 0; sbottom = 0; //'sun//'s down
                                                                                    //'night = "y" //'it//'s officially night
      } //Endif//
        //'now compute long wave irradiance
      ea = patm * (maxvpd - vpd); //'vapor pressure in kPa
                                  //'ta = Cells(3, 19) + 273.15 //'air temp in K
      eac = 1.72 * pow((ea / (airtemp + 273.15)), (1.0 / 7.0)); //'emissivity of clear sky CN  10.10
                                                                //'boltz = Cells(9, 11) //'boltzman constant
      la = eac * sbc * pow((airtemp + 273.15), 4); //'long wave irradiance from clear sky
      lg = 0.97 * sbc * pow((airtemp + 273.15), 4); //'long wave irradiance from ground...assumes equilibrium with air temp
                                                    //'Cells(14 + i, 16) = la
                                                    //'Cells(14 + i, 17) = lg
   }

   void getsoilwetness() //'gets predawns accounting for rain, groundwater flow, transpiration, and redistribution via soil and roots
   {
      double tempDouble = 0.0;

      drainage = 0;
      runoff = 0;
      if (dd == 1 || isNewYear) { //if// //'every layer starts at initial % of field capacity

         if (!(useGSData && gs_yearIndex > 0))
         {
            waterold = 0; //'total root zone water content (mmol m-2 ground area)
            for (z = 0; z <= layers; z++)//z = 0 To layers //'
            {
               x = 10; //'MPa water potential for getting residual thetafrac
               thetafracres[z] = swc(x); //'residual thetafrac
               thetafracfc[z] = (1 - thetafracres[z]) * fieldcapfrac + thetafracres[z]; //'thetafrac at field capacity
               thetafc[z] = thetasat[z] * thetafracfc[z]; //'water content at field capacity
            } //for//z //'
            for (z = 0; z <= layers; z++)//z = 0 To layers
            {
               water[z] = thetafc[z] * depth[z]; //'field capacity estimated as 1/2 saturated capacity, water content of layer in m3 water per m2 ground area
               fc[z] = water[z]; //'records field capacity in m3 water volume per m2 ground area.
               water[z] = ffc * water[z]; //'start off with initial fraction of field capacity
                                          //'Cells(16 + dd, 42 + z) = water[z]
               waterold = waterold + water[z]; //'in m3/m2 ground area
            } //for//z
            dSheet.Cells(rowD + dd, colD + dColF_End_watercontent) = waterold * 1000; //'root zone water content in mm m-2 ground area
                                                                                      //'waterold = waterold * 55555556# //'convert m3 water per ground area to mmol water per ground area
                                                                                      // [HNT] starting water now counts as an input
            gs_ar_input[gs_yearIndex] = gs_ar_input[gs_yearIndex] + waterold * 1000;
         }

         if (gs_yearIndex == 0) // if it's the first year, store this as the off-season starting water because we don't have a real value
            gs_ar_waterInitial_OFF[gs_yearIndex] = waterold * 1000;

         // store the initial water content to check how much we consume at the end
         gs_ar_waterInitial[gs_yearIndex] = waterold * 1000;
         // [/HNT]
      } //End if// //'dd=1 if
        //'if pet = "y" Or pet = "n" { //if// //'do the original routine...doesn//'t run for PET scenario
      if ((dd > 1 && !isNewYear) || (useGSData && gs_yearIndex > 0)) { //if// //'get flows that happened during previous timestep
         for (z = 0; z < layers; z++)//z = 0 To layers - 1 //'transpiration, root and soil redistribution
         {
            if (night == "n") { //if// //'it//'s day, must adjust elayer for sun vs. shade weighting
               layerflow = elayer[z][halt] * laisl / lai + elayer[z][haltsh] * laish / lai; //'weighted flow; NOTE: elayer = 0 for layer 0
            }
            else {
               layerflow = elayer[z][halt]; //'no adjustment necessary at night; NOTE: elayer = 0 for layer 0
            } //End if// //'night if
            layerflow = layerflow * baperga * 1.0 / 998.2 * timestep; //'rootflow into (= negative rootflow) or out (positive flow) of layer in m3/m2 ground area
            layerflow = layerflow + soilredist[z] * 1.0 / 998.2 * timestep; //'redistribution between layers (negative is inflow, positive is outflow). NOTE: soilredist(0) includes soil evaporation for layer 0
            water[z] = water[z] - layerflow; //'subtracts rootflow from layer on per ground area basis
         } //for//z
           //'now do the bottom layer and potential groundwater input
         if (night == "n") { //if// //'it//'s day, must adjust layerflow for sun vs. shade weighting
            layerflow = elayer[layers][halt] * laisl / lai + elayer[layers][haltsh] * laish / lai; //'weighted flow
         }
         else {
            layerflow = elayer[layers][halt]; //'no adjustment necessary at night
         } //End if// //'night if
         layerflow = layerflow * baperga * 1 / 998.2 * timestep; //'rootflow into (= negative rootflow) or out (positive flow) of layer in m3/m2 ground area
         layerflow = layerflow + soilredist[layers] * 1 / 998.2 * timestep; //'redistribution between layers (negative is inflow, positive is outflow)
                                                                            //'water(layers) = water(layers) - layerflow //'subtracts rootflow from layer on per ground area basis
         if (layerflow < 0) { //if// //'water is added
            for (z = layers; z >= 0; z--)//z = layers To 0 Step -1 //'start at bottom, go up
            {
               deficit = thetasat[z] * depth[z] - water[z]; //'m of water required to wet up layer to SATURATION
               if (-1 * layerflow - deficit >= 0) { //if// //'there//'s enough to wet the layer...remember, negative flow is flow into the layer
                  water[z] = thetasat[z] * depth[z]; //'m water at saturation in layer
                  layerflow = layerflow + deficit; //'reduce what//'s left over for next layers
               }
               else { //'just soak up all the groundwater
                  water[z] = water[z] - layerflow; //'add to bottom layer
                  layerflow = 0; //' all gone
               } //End if// //'wetting if
            } //for//z
            runoff = runoff - layerflow; //'add what//'s left to runoff...
         }
         else { //'groundwater is positive...bottom layer is losing water
            water[layers] = water[layers] - layerflow; //'subtract from bottom layer
         } //End if// //'layerflow if

           //'now reset any exhausted layers to extraction limit
         if (water[0] <= 0)
            water[0] = 0.00001; //'set lower limit to surface water content
         for (z = 1; z <= layers; z++)//z = 1 To layers
         {
            if (water[z] < swclimit[z])
               water[z] = swclimit[z]; //'water at limit
         } //for//z

         bool rainOverride = false;
         if (iter_Counter == 0 && tod == 23 && (stage_ID == STAGE_ID_HIST_OPT || stage_ID == STAGE_ID_FUT_OPT))
            rainOverride = true;

         //'now check for rain during PREVIOUS TIME STEP
         if (dSheet.Cells(rowD + dd - 1, colD + dColRain) >= 0.0 || rainOverride == true) { //if// //'there//'s been rain!
            rain = dSheet.Cells(rowD + dd - 1, colD + dColRain) * 0.001; //'convert mm of rain to m depth...equivalent to m3/m2 water volume per ground area
                                                                         //'if raining(xx) = "n" { //if// rain = 0 //'rain turned off

            if (rainOverride == true)
               rain = 0.1; // 100 mm of rain every night if we're doing a BAGA optimization reference K run (iter_Counter is zero)

            if (rainEnabled == false)
               rain = 0;
            // [HNT] temp
            //if (rain > 0.1) // more than 1/10th meter per hour is a data anomoly, however this is a poor assumption and these should be checked in the weather curator
            //   rain = 0.0;
            // [/HNT]
            //'sumrain = rain * 1000 + sumrain //'total rain input in mm m-2
            for (z = 0; z <= layers; z++)//z = 0 To layers //'add in rain
            {
               if (rain <= 0)
                  break; //'rain//'s used up
               deficit = fc[z] - water[z]; //'m3 of water per m2 to wet up layer to fc
                                           //if (deficit >= -0.001) { //if// //'layer is at or below field capacity
               if (rain - deficit >= 0) { //if// //'there//'s enough to wet the layer
                  water[z] = fc[z];
                  rain = rain - deficit;
                  //'sumrain = sumrain + deficit //'absorbed rain
                  drainage = rain * 1000.0; //'any left over will drain out the bottom unless layer is rising above field capacity
               }
               else {
                  water[z] = water[z] + rain; //'it all goes into the first layer
                                              //'sumrain = sumrain + rain
                  rain = 0; //'rain used up
                  drainage = 0;
               } //End if// //'wetting up to field capacity "if"
                 //}
            } //for//z

              // If there's drainage, and the ground water is on, that should be used to fill up to saturation
            if (rain > 0.0 && ground == "y") // the remaining "drainage" is also still stored in the rain variable
            {
               // this is kind of inefficient, but the rain routine actually drained all the layers to FC even if GW was on and we should have been filling to sat
               // now we start at the bottom and fill the layers to saturation using the drainage

               for (j = layers; j >= 0; j--) //j = z - 1 To 0 Step -1 //'go back up to fill profile to saturation
               {
                  if (rain <= 0) //if// Exit for
                     break;
                  deficit = thetasat[j] * depth[j] - water[j];
                  if (deficit >= 0) { //if// //'got capacity
                     if (rain - deficit >= 0) { //if// //'enough rain to saturate the layer
                        water[j] = thetasat[j] * depth[j]; //'saturate the layer
                        rain = rain - deficit; //'reduce rain
                     }
                     else { //'rain absorbed by layer
                        water[j] = water[j] + rain;
                        rain = 0; //'use up rain
                     } //End if// //'deficit=>0 "if"
                  }
                  else { //'deficit<0...layer//'s saturated
                     rain = rain - deficit; //'increase rain by super-saturated amount (deficit is negative)
                     water[j] = thetasat[j] * depth[j]; //'reset to saturation
                  } //End if// //'deficit <>0 if
               } //for//j
               runoff = runoff + rain; //'whatever is left over will run off
               drainage = 0; //'no drainage if any layer is rising above field capacity
               rain = 0; //'reset rain to zero
            }
            //'sumdrain = sumdrain + drainage //'total drainage
         } //End if// //'rain if

           //'now check for exhausted layers
         if (water[0] <= 0)
            water[0] = 0.00001; //'set lower limit to surface water content
         for (z = 1; z <= layers; z++)//z = 1 To layers
         {
            if (water[z] < swclimit[z])
               layer[z] = 1; //'water exhausted
         } //for//z
           //'now get water content change over PREVIOUS time step
         if ((dd > 1 && !isNewYear) || (useGSData && gs_yearIndex > 0)) { //if// //'now get updated new water content
            waternew = 0;
            for (z = 0; z <= layers; z++)//z = 0 To layers //'check for exhausted layers
            {
               waternew = water[z] + waternew;
               //'Cells(16 + dd, 42 + z) = water[z]
            } //for//z



              //'waternew = waternew * 55555556# //'new water content at beginning of timestep
            waterchange = waternew - waterold; //'total water change in mm m-2 ground area
                                               //'waterchange = waterchange * 1 / baperga * 1 / laperba //'total water in mmol per m2 leaf area
            dSheet.Cells(rowD + dd, colD + dColF_End_watercontent) = waternew * 1000; //'root zone water content in mm

                                                                                      // always store the water content as the "final" -- not worth testing if it's really the last day of the year
            gs_ar_waterFinal[gs_yearIndex] = waternew * 1000;

            dSheet.Cells(rowD + dd, colD + dColF_End_waterchange) = waterchange * 1000; //'change in water content over PREVIOUS timestep
                                                                                        //'if raining(xx) = "y" { //if// Cells(16 + dd, 59) = Cells(16 + dd - 1, 4) //'rain input per previous timestep
            if (rainEnabled == true)// && dSheet.Cells(rowD + dd - 1, colD + dColRain) < 100.0)
               dSheet.Cells(rowD + dd, colD + dColF_End_rain) = dSheet.Cells(rowD + dd - 1, colD + dColRain); //'rain input per previous timestep
            dSheet.Cells(rowD + dd, colD + dColF_End_gwater) = gwflow; //'groundwater input in mm per timestep
                                                                       //'Cells(16 + dd, 61) = transpirationtree * 3600 * timestep * laperba * baperga * 0.000000018 * 1000 //'transpiration per ground area in mm m-2 per timestop
            dSheet.Cells(rowD + dd, colD + dColF_End_drainage) = drainage; //'total drainage in mm per timestep
            dSheet.Cells(rowD + dd, colD + dColF_End_input) = dSheet.Cells(rowD + dd, colD + dColF_End_rain) + dSheet.Cells(rowD + dd, colD + dColF_End_gwater); //'total input per timestep in mm
                                                                                                                                                                 //'Cells(16 + dd, 64) = water(0) * 1000 //'water in top layer in mm
            gs_ar_input[gs_yearIndex] = gs_ar_input[gs_yearIndex] + dSheet.Cells(rowD + dd, colD + dColF_End_input);

            dSheet.Cells(rowD + dd, colD + dColF_End_runoff) = runoff * 1000; //'excess root zone water per timestep in mm
            waterold = waternew; //'reset to beginning of current timestep
         } //End if// //'dd>1 if
      } //End if// //'dd>1 if
        //'} //End if////'pet if
      if (dd > 1 && !isNewYear) { //if//
         tempDouble = transpirationtree * 3600 * timestep * laperba * baperga * 0.000000018 * 1000;
         dSheet.Cells(rowD + dd, o + colD + dColF_End_E) = tempDouble;//transpirationtree * 3600 * timestep * laperba * baperga * 0.000000018 * 1000; //'transpiration per ground area in mm m-2 per timestop
         if (gs_inGrowSeason) // only record growing season E (should not be any non-GS E, but just for safety)
            gs_ar_E[gs_yearIndex] = gs_ar_E[gs_yearIndex] + tempDouble;

         dSheet.Cells(rowD + dd, o + colD + dColF_End_soilEvap) = soilevap * 1 / 998.2 * timestep * 1000; //'evaporative water loss in mm per timestep

         tempDouble = transpirationtree * 3600 * timestep * laperba * baperga * 0.000000018 * 1000 + soilevap * 1 / 998.2 * timestep * 1000;
         dSheet.Cells(rowD + dd, o + colD + dColF_End_ET) = tempDouble;//transpirationtree * 3600 * timestep * laperba * baperga * 0.000000018 * 1000 + soilevap * 1 / 998.2 * timestep * 1000;
         gs_ar_ET[gs_yearIndex] = gs_ar_ET[gs_yearIndex] + tempDouble;

         tempDouble = atree * timestep * 3600 * 0.001;
         dSheet.Cells(rowD + dd, o + colD + dColF_End_ANet) = tempDouble;//atree * timestep * 3600 * 0.001; //'Anet per timestep in mmoles per leaf area
         if (!std::isnan(tempDouble) && gs_inGrowSeason) // only record growing season A
         {
            gs_ar_Anet[gs_yearIndex] = gs_ar_Anet[gs_yearIndex] + tempDouble;
            //anytime we record A, also record the ci-related outputs
            // TODO CRIT units???
            // Anet in mmoles per leaf area as in final output columns? (calc above)
            if (night == "n") // daytime only
            {
               gs_ar_cica[gs_yearIndex] += cinc / ca; // these are in mols/mol, but it's a ratio so not important
               gs_ar_cica_N[gs_yearIndex]++;

               gs_ar_Aci[gs_yearIndex] += tempDouble * cinc; // for Aci what units? pa would be (cinc * patm * 1000)
               gs_ar_AnetDay[gs_yearIndex] += tempDouble; // keep a daytime-only Anet tally
                                                          //this is for A-weighted Ci
                                                          //gs_ar_Acica[gs_yearIndex] += tempDouble * cinc / ca; // for this one, another ratio so ignore units
            }
         }
      } //End if// //'dd>1 if

      if (tod == 16 && !gs_doneFirstDay && gs_inGrowSeason && dSheet.Cells(rowD + dd - 3, colD + dColF_CP_kplant) > 0.000000001) { //if// //'get midday k//'s for day 1
                                                                                                                                   // VPD zero case -- if the stomata did not open on the first day of the GS, kplant won't have been set and will be zero... in which case, try again tomorrow
         gs_doneFirstDay = true;

         sum = 0;
         for (z = 1; z <= 3; z++)//z = 1 To 3
         {
            sum = sum + dSheet.Cells(rowD + dd - z, colD + dColF_CP_kplant);
         } //for//z
         kpday1 = sum / 3.0; //'average midday kplant on day one
         sum = 0;
         for (z = 1; z <= 3; z++)//z = 1 To 3
         {
            sum = sum + dSheet.Cells(rowD + dd - z, colD + dColF_CP_kxylem);
         } //for//z
         kxday1 = sum / 3.0; //'average midday kxylem on day one

                             // [HNT] first day of the new growing season! Record the starting water
         gs_ar_waterInitial_GS[gs_yearIndex] = dSheet.Cells(rowD + dd, colD + dColF_End_watercontent); // no matter what happened above, this has been set by this point
                                                                                                       // and record this as the FINAL water for THIS YEAR's off season -- remember that this year's off season is the PRECEDING winter
         gs_ar_waterFinal_OFF[gs_yearIndex] = dSheet.Cells(rowD + dd, colD + dColF_End_watercontent);
         // [/HNT]
      } //End if// //'dd=16 if
      if (gs_doneFirstDay) { //if// //'calculate plc relative to midday of day 1
         if (iter_refK < 0.000000001) // || iter_Counter == 0 // no longer appropriate to test for iter_Counter == 0 ... may be doing a seperate stress profile that refers to a saved refK, which will have been loaded at start of modelProgramMain
            tempDouble = 100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kplant) / kpday1); //' done to avoid repeating this calculation
         else // if we haven't loaded a ref K it will be set to zero, and we need to fall back to the old method.. otherwise use refK to calculate PLC
         {
            tempDouble = 100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kplant) / iter_refK); // PLCp calculated from refK if it exists
            if (tempDouble < 0.0) // if we're using the refK, it's possible for this to be negative briefly -- should be considered zero
               tempDouble = 0.0;
         }
         dSheet.Cells(rowD + dd, colD + dColF_End_PLCplant) = tempDouble; //'100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kplant) / kpday1) 'plc plant...prior timestep

                                                                          // no matter what, add it to the tally for calculating the mean over-season PLCp
         gs_ar_PLCSum[gs_yearIndex] = gs_ar_PLCSum[gs_yearIndex] + tempDouble; // yearly running total of PLC values
         gs_ar_PLCSum_N[gs_yearIndex] = gs_ar_PLCSum_N[gs_yearIndex] + 1; // total hours in GS
                                                                          // now test for highest PLC and hours > 85
         if (tempDouble > gs_ar_PLCp[gs_yearIndex])
            gs_ar_PLCp[gs_yearIndex] = tempDouble;
         if (tempDouble > 85.0)
            gs_ar_PLC85[gs_yearIndex] = gs_ar_PLC85[gs_yearIndex] + 1;

         tempDouble = 100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kxylem) / kxday1);
         dSheet.Cells(rowD + dd, colD + dColF_End_PLCxylem) = tempDouble; //'100 * (1 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kxylem) / kxday1) 'plc xylem...prior timestep
         if (tempDouble > gs_ar_PLCx[gs_yearIndex])
            gs_ar_PLCx[gs_yearIndex] = tempDouble;

         //dSheet.Cells(rowD + dd, colD + dColF_End_PLCplant) = (100.0 * (1.0 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kplant) / kpday1)); //'plc plant...prior timestep
         //dSheet.Cells(rowD + dd, colD + dColF_End_PLCxylem) = (100.0 * (1.0 - dSheet.Cells(rowD - 1 + dd, colD + dColF_CP_kxylem) / kxday1)); //'plc xylem...prior timestep

         // [HNT] keep track of in-season input
         if (gs_inGrowSeason)
         {
            // if we've done the first GS day and we're in the growing season then it's "this year" during GS
            gs_ar_waterInput_GS[gs_yearIndex] += dSheet.Cells(rowD + dd, colD + dColF_End_input);
         }
         else
         {
            // if we've done the first GS day but we're NOT in the growing season, it's the winter following GS.
            // this is considered NEXT YEAR's off-season input! First check if we already have a value for the FINAL water for THIS YEAR, because if it's zero then
            // this is the first timestep of the winter and we need to store it
            if (gs_ar_waterFinal_GS[gs_yearIndex] <= 0.0 && gs_ar_waterInitial_OFF[gs_yearIndex + 1] <= 0.0) // ok to use == or <= with double here because this will be memset to 0 if it hasn't been set
            {
               gs_ar_waterFinal_GS[gs_yearIndex] = dSheet.Cells(rowD + dd, colD + dColF_End_watercontent); // ok to overshoot by 1 hour, I think (instead of using dd - 1)
               gs_ar_waterInitial_OFF[gs_yearIndex + 1] = dSheet.Cells(rowD + dd, colD + dColF_End_watercontent); // also the off-season initial for next year
            }
            else // otherwise, we're in the middle of NEXT YEAR's off season ... note this +1 on an array index is super lazy and bad. Make sure to never run exactly this array size # of years
            {
               gs_ar_waterInput_OFF[gs_yearIndex + 1] += dSheet.Cells(rowD + dd, colD + dColF_End_input); // add the stored input to the input tally, for NEXT YEAR
            }
         }
         // [/HNT]
      } //End if// //'dd>16 if
      else if (!gs_inGrowSeason)// have NOT done first day and are NOT in growing season
      {
         // we must be in the pre-GS winter of what we called the NEXT year above... so gs_yearIndex is now that year, and we add to the off-season input for THIS year
         gs_ar_waterInput_OFF[gs_yearIndex] += dSheet.Cells(rowD + dd, colD + dColF_End_input);
      }
   }

   void getpredawns()
   {
      //'first check for layer participation...only rooted layers, not layer 0

      for (k = 1; k <= layers; k++)//k = 1 To layers //'assign source pressures, set layer participation
      {
         if (layerfailure[k] == "root") { //if//
            if (kminroot[k] == 0)
               layer[k] = 1; //'gone from scene if roots cavitated at midday
            if (kminroot[k] != 0) { //if// //'roots still around
               layerfailure[k] = "ok";
               layer[k] = 0;
            } //End if//
         } //End if//
         if (layerfailure[k] == "rhizosphere")
            layer[k] = 0; //'layer can come back to life
      } //for//k
        //'after getting water[z] and layer participation, get predawns
      for (z = 0; z <= layers; z++)//z = 0 To layers
      {
         if (layer[z] == 0) { //if//
            theta = water[z] / depth[z]; //'convert m3 water per m2 ground back to m3 water / m3 soil
            x = theta / thetasat[z]; //'remember, VG function takes theta/thetasat as input
            pd[z] = rvg(x); //'soil pressure of layer
            prh[z] = pd[z]; //'guess for NR solution
            if (pd[z] >= pcritrh[z] && z > 0) { //if// //'only rooted layers // [HNT] >= instead of > for consistency w/ Newton Rhapson update
               layer[z] = 1;
               layerfailure[z] = "rhizosphere";
            } //End if//
            if (pd[z] >= pcritr[z] && z > 0) { //if// //'only rooted layers // [HNT] >= instead of > for consistency w/ Newton Rhapson update
               layer[z] = 1;
               layerfailure[z] = "root";
               kminroot[z] = 0;
            } //End if//
         }
         else { //'layer//'s disconnected
            pd[z] = pcritr[z];
            prh[z] = pcritr[z];
         } //End if//
      } //for//z
        //'now get guess of proot
      sum = 0;
      t = 0;
      for (k = 1; k <= layers; k++)//k = 1 To layers
      {
         if (layer[k] == 0) { //if//
            sum = sum + pd[k];
         }
         else { //'predawn is not seen by the roots
            t = t + 1;
         } //End if//
      } //for//k
      failspot = "no failure";
      if (t < layers) { //if//
         pr = sum / (layers - t); //'set unknown proot to average pd
         prinitial = pr; //'store initial value if NR gets off the rails
      }
      else {
         failure = 1;
      } //End if//
      for (z = 1; z <= layers; z++)//z = 1 To layers
      {
         dSheet.Cells(rowD + dd, colD + dColF_p1 - 1 + o + z) = pd[z]; //'soil pressures by layer (only for rooted layers)
      } //for//z
        //'Cells(16 + dd, 65) = pd(0) //'water potential of top layer
   }

   void soilflow() //'gets flow out of each layer via soil, not including the groundwater basement
   {
      for (z = 0; z < layers; z++)//z = 0 To layers - 1
      {
         store = kmaxrh[z]; //'store the rhizosphere kmax for now
         kmaxrh[z] = kkmax[z + 1] * 1 / (depth[z] / 2.0 + depth[z + 1] / 2.0); //'reset it to the vertical soil kmax per m2 ground using distance between midpoints of each layer. Use properties of fatter layer.
         if (pd[z] == pd[z + 1])
            soilredist[z] = 0; //'no redistribution (gravity ignored!)
         if (pd[z] != pd[z + 1]) { //if//
            if (pd[z] < pd[z + 1]) { //if// //'flow is out of layer z (positive)
               p1 = pd[z];
               pend = pd[z + 1];
            }
            else { //'flow is into layer z(negative)
               p1 = pd[z + 1];
               pend = pd[z];
            } //End if//
            e = 0; //'flow integral
            do
            {
               p2 = p1 + pinc;
               qtrapvg(p1, p2, s);
               e = e + s;
               p1 = p2; //'reset p1 for next increment
            } while (!(p2 >= pend || p2 > 5));
            //Loop Until p2 >= pend Or p2 > 5 //'cutoff for really dry soil
            if (pd[z] < pd[z + 1]) { //if// //'flow is out of layer z (positive)
               soilf[z] = e; //'flow in kg hr-1 m-2 ground area
            }
            else { //'flow is into layer z(negative)
               soilf[z] = -e;
            } //End if//
         } //End if//
         kmaxrh[z] = store; //'reset rhizosphere kmaxrh (though i don//'t THINK it//'s used again?)
      } //for//z
      soilredist[0] = soilf[0]; //'set upper layer which has only one flux
      soilredist[layers] = -1 * soilf[layers - 1]; //'set water flowing into/out of top of the bottom layer...soilredist is net soil flow of layer
                                                   //'now calculate net flows of internal layers
      for (z = 1; z < layers; z++)//z = 1 To layers - 1
      {
         soilredist[z] = -1 * soilf[z - 1] + soilf[z]; //'add up to get net flow...remember negative flow is flow into the layer
      } //for//z
   }

   void deepflow() //'gets flow into/out of the lowermost layer
   {
      z = layers; //'just the bottom layer
      store = kmaxrh[z]; //'store the rhizosphere kmax for now
      kmaxrh[z] = kkmax[z] / grounddistance; //'reset it to the vertical soil kmax using distance to groundwater
      if (pd[z] >= pground) { //if// //'groundwater will flow in (negative flow)
         p1 = pground;
         pend = pd[z];
         e = 0; //'flow integral
         do
         {
            p2 = p1 + pinc;
            qtrapvg(p1, p2, s);
            e = e + s;
            p1 = p2; //'reset p1 for next increment
         } while (!(p2 >= pend));
         //Loop Until p2 >= pend
         groundflow = -1 * e; //'negative flow (into bottom layer)
         gwflow = -1 * groundflow * 1 / 998.2 * timestep * 1000; //'groundwater flow totals in mm m-2 ground area per time step
      }
      else { //'flow is out of bottom layer(positive)
         p1 = pd[z];
         pend = pground;
         e = 0;//'flow integral
         do
         {
            p2 = p1 + pinc;
            qtrapvg(p1, p2, s);
            e = e + s;
            p1 = p2; //'reset p1 for next increment
         } while (!(p2 >= pend));
         //Loop Until p2 >= pend
         groundflow = e; //'positive flow (out of bottom layer)
         drainage = drainage + groundflow * 1 / 998.2 * timestep * 1000; //'drainage totals in mm m-2 ground area per time step
      } //End if//
      kmaxrh[z] = store; //'reset rhizosphere kmaxrh (though i don//'t THINK it//'s used again?)
                         //'now calculate net flow into lowermost layer
      soilredist[z] = soilredist[z] + groundflow; //'add up to get net flow...remember negative flow is flow into the layer
   }

   void soilevaporation() //'get soil evaporation
   {
      emission = emiss * sbc * pow((soiltemp + 273.15), 4); //'long wave emissivity of soil
      lc = 0.97 * sbc * pow((leaftempsh[haltsh] + 273.15), 4); //'thermal emissivity from bottom of canopy
      rabssoil = soilabssol * sbottom + 0.97 * lc; //'rabs of soil surface; does not account for slope (beam only)
      mdensair = 44.6 * (patm / 101.3) * 273.15 / (airtemp + 273.15); //'molar density of air, CN 3.3, moles/m3
      gha = 0.16 * mdensair * us / (log((xh - zdispl) / rough) * log((xh - zdispl) / zh)); //'CN 7.28/14.9,heat and vapor aerodynamic conductance for soil, mol m-2 s-1, stability factors = 0
      soilep = (rabssoil - emission - sha * (soiltemp - airtemp) * (grad + gha)) / lambda; //'potential evaporation in moles m-2 s-1; 12.8
      if (soilep < 0)
         soilep = 0; //'set to zero if negative
      rha = 1 - vpd / maxvpd; //'relative humidity of air
      rhs = exp((-pd[0] * 0.00001805) / (gas * 0.000001 * (soiltemp + 273.15))); //'relative humidity in equilibrium with soil surface water potential
      if (rha == 1) { //if//
         soilevap = 0;
      }
      else {
         soilevap = soilep * (rhs - rha) / (1 - rha); //'campbell 1985, eqn 9.14; actual soil evaporation rate in moles m-2s-1
      } //End if//
      if (soilevap < 0)
         soilevap = 0; //'don//'t let it go negative
      soilevap = soilevap * (1 - baperga); //'reduce for basal area
      soilevap = soilevap * 0.0180153 * 3600; //'converts from moles m-2s-1 to kg m-2 hr-1
      soilredist[0] = soilredist[0] + soilevap; //'add evaporative loss to redistribution for top layer
   }

   void leaftemps()  //'gets sun layer leaf temp and leaf-to-air vpd from E, Campbell and Norman book
   {
      rabs = 0.5 * (0.5 * ssun + 0.5 * sref) + emiss * (0.5 * la + 0.5 * lg); //'total absorbed radiation for sun leaves; CN 11.14
      lambda = -42.9143 * airtemp + 45064.3; //'heat of vaporization for water at air temp in J mol-1
      grad = 0.1579 + 0.0017 * airtemp + 0.00000717 * pow(airtemp, 2); //'radiative conductance (long wave) at air temp in mol m-2 s-1
      gha = 1.4 * 0.135 * pow((wind / leafwidth), 0.5); //'heat conductance in mol m-2s-1
      eplantl[p] = eplant[p] * (1.0 / laperba) * (1.0 / 3600.0) * 55.4; //'convert to E per leaf area in mol m-2s-1
      numerator = rabs - emiss * sbc * pow((airtemp + 273.2), 4) - lambda * eplantl[p] / 2.0; //'divide E by 2 because energy balance is two sided.
      denominator = sha * (grad + gha);
      leaftemp[p] = airtemp + numerator / denominator; //'leaf temp for supply function
      lavpd[p] = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
      lavpd[p] = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * leaftemp[p])) - lavpd[p] + vpd; //'leaf-to-air vpd
      if (lavpd[p] < 0)
         lavpd[p] = 0; //'don//'t allow negative lavpd
                       //'eplantl[p] = eplantl[p] * 1000 //'convert to mmol m-2 s-1
   }

   void leaftempsshade()  //'gets shade layer leaf temp and leaf-to-air vpd from E, Campbell and Norman book
   {
      rabs = 0.5 * (0.5 * sshade + 0.5 * sref) + emiss * lg; //'total absorbed radiation for shaded leaves
                                                             //'lambda = -42.9143 * airtemp + 45064.3 //'heat of vaporization for water at air temp in J mol-1
                                                             //'grad = 0.1579 + 0.0017 * airtemp + 0.00000717 * airtemp ^ 2 //'radiative conductance (long wave) at air temp in mol m-2 s-1
                                                             //'gha = 1.4 * 0.135 * (wind / leafwidth) ^ 0.5 //'heat conductance in mol m-2s-1
                                                             //'eplantl[p] = eplant[p] * (1 / laperba) * (1 / 3600) * 55.4 //'convert to E per leaf area in mol m-2s-1
      numerator = rabs - emiss * sbc * pow((airtemp + 273.2), 4) - lambda * eplantl[p] / 2.0; //'divide E by 2 because energy balance is two sided.
      denominator = sha * (grad + gha);
      leaftempsh[p] = airtemp + numerator / denominator; //'leaf temp for supply function
      lavpdsh[p] = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
      lavpdsh[p] = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * leaftempsh[p])) - lavpdsh[p] + vpd; //'leaf-to-air vpd
      if (lavpdsh[p] < 0)
         lavpdsh[p] = 0; //'don//'t allow negative lavpd
      eplantl[p] = eplantl[p] * 1000; //'convert to mmol m-2 s-1
   }

   void assimilation() //'gets assimilation for sun leaves
   {
      //'get g from D and e
      if (lavpd[p] == 0) { //if//
         gcanw[p] = gmax; //'maximum g if no lavpd
      }
      else {
         gcanw[p] = eplantl[p] / lavpd[p]; //'gcanopy in mmol m-2s-1 (leaf area)
      } //End if//
        //'gcanw[p] = eplantl[p] / lavpd[p] //'gcanopy in mmol m-2s-1 (leaf area)
      gcanc[p] = (gcanw[p] / 1.6) * 1000; //'convert to CO2 conductance in umol m-2 s-1
                                          //'adjust photosynthetic inputs for Tleaf
      comp = comp25 * exp((37830 * ((leaftemp[p] + 273.15) - 298.15)) / (298.15 * gas * (leaftemp[p] + 273.15))); //'Bernacchi via Medlyn
      numerator = (1 + exp((svvmax * 298.2 - hdvmax) / (gas * 298.2))) * exp((havmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftemp[p])));
      denominator = 1 + exp((svvmax * (273.2 + leaftemp[p]) - hdvmax) / (gas * (273.2 + leaftemp[p])));
      vmax = vmax25 * numerator / denominator; //'vmax corrected via Leunig 2002
      numerator = (1 + exp((svjmax * 298.2 - hdjmax) / (gas * 298.2))) * exp((hajmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftemp[p])));
      denominator = 1 + exp((svjmax * (273.2 + leaftemp[p]) - hdjmax) / (gas * (273.2 + leaftemp[p])));
      jmax = jmax25 * numerator / denominator; //'jmax corrected via Leunig 2002
      kc = kc25 * exp((79430 * ((leaftemp[p] + 273.15) - 298.15)) / (298.15 * gas * (leaftemp[p] + 273.15))); //'Bernacchi via Medlyn
      ko = ko25 * exp((36380 * ((leaftemp[p] + 273.15) - 298.15)) / (298.15 * gas * (leaftemp[p] + 273.15))); //'Bernacchi via Medlyn
      rday25 = vmax25 * 0.01; //'from Medlyn 2002
      rday[p] = rday25 * pow(2, ((leaftemp[p] - 25) / 10.0));
      rday[p] = rday[p] * pow((1 + exp(1.3 * (leaftemp[p] - 55))), -1); //'high temp inhibition, collatz
      if (night == "n" && gcanc[p] > 0) { //if// //'solve for A and ci
         if (p == 1) { //if// //'stomata have just opened: find the mitochondrial compensation point
            ci = comp - 0.00000001; //'start at photorespiratory compensation point
            do //'loop through A-ci curve to find mitochondrial compensation point
            {
               ci = ci + 0.00000001;
               jact = ((qmax * qsl + jmax) - pow((pow((-qmax * qsl - jmax), 2) - 4 * lightcurv * qmax * qsl * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
               je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
               numerator = vmax * (ci - comp);
               denominator = ci + kc * (1 + oa / ko);
               jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
               var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
               var = var - rday[p]; //'convert gross to NET A
            } while (!(var >= 0 || ci >= ca));
            //Loop Until var >= 0 Or ci >= ca //'there//'s positive Anet or not enough light
            psyn[p] = var; //'always start with zero or above
            cin[p] = ci; //'the dark compensation point...not predicting negative Anet
         } //End if// //'p=1 if
         if (p > 1) { //if//
            ci = cin[p - 1]; //'start from previous ci
            do //'loop through A-ci curve to find ci
            {
               ci = ci + 0.00000001;
               marker = gcanc[p] * (ca - ci); //'marker is NET A from g (in umol m-2s-1, need to convert ca-ci to mole fraction), gets smaller as ci increase
               jact = ((qmax * qsl + jmax) - pow((pow((-qmax * qsl - jmax), 2) - 4 * lightcurv * qmax * qsl * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
               je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
               numerator = vmax * (ci - comp);
               denominator = ci + kc * (1 + oa / ko);
               jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
               var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
               var = var - rday[p]; //'convert gross to NET A
            } while (!(var >= marker || ci >= ca));
            //Loop Until var >= marker Or ci >= ca //'when "a-ci" value just reaches the "g" value, that//'s the right ci...
            psyn[p] = var; //'net assimilation, umol m-2 leaf area s-1
            if (psyn[p] > psynmax)
               psynmax = psyn[p];
            cin[p] = ci; //'store ci
         } //End if// //'p>1 if
      }
      else { //'it//'s night or stomata are closed
         if (night == "y")
         {
            psyn[p] = 0 - rday[p];
            psynmax = 0;
            cin[p] = ca; //'respiration accounted for at night
         }
         if (night == "n") {
            psyn[p] = 0;
            psynmax = 0;
            cin[p] = ca; //'it//'s day and p=0, g=0
         }
      } //End if// //'night if
        //'Cells(16 + p, 65) = psynsh[p]
   }

   void assimilationshade() //'gets assimilation for shade leaves
   {
      //'get g from D and e
      if (lavpdsh[p] == 0) { //if//
         gcanwsh[p] = gmax; //'set to maximum if no vpd
      }
      else {
         gcanwsh[p] = eplantl[p] / lavpdsh[p]; //'gcanopy in mmol m-2s-1 (leaf area)
      } //End if//
        //'gcanwsh[p] = eplantl[p] / lavpdsh[p] //'gcanopy in mmol m-2s-1 (leaf area)
      gcancsh[p] = (gcanwsh[p] / 1.6) * 1000; //'convert to CO2 conductance in umol m-2 s-1
                                              //'adjust photosynthetic inputs for Tleaf
      comp = comp25 * exp((37830 * ((leaftempsh[p] + 273.15) - 298.15)) / (298.15 * gas * (leaftempsh[p] + 273.15))); //'Bernacchi via Medlyn
      numerator = (1 + exp((svvmax * 298.2 - hdvmax) / (gas * 298.2))) * exp((havmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftempsh[p])));
      denominator = 1 + exp((svvmax * (273.2 + leaftempsh[p]) - hdvmax) / (gas * (273.2 + leaftempsh[p])));
      vmax = vmax25 * numerator / denominator; //'vmax corrected via Leunig 2002
      numerator = (1 + exp((svjmax * 298.2 - hdjmax) / (gas * 298.2))) * exp((hajmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftempsh[p])));
      denominator = 1 + exp((svjmax * (273.2 + leaftempsh[p]) - hdjmax) / (gas * (273.2 + leaftempsh[p])));
      jmax = jmax25 * numerator / denominator; //'jmax corrected via Leunig 2002
      kc = kc25 * exp((79430 * ((leaftempsh[p] + 273.15) - 298.15)) / (298.15 * gas * (leaftempsh[p] + 273.15))); //'Bernacchi via Medlyn
      ko = ko25 * exp((36380 * ((leaftempsh[p] + 273.15) - 298.15)) / (298.15 * gas * (leaftempsh[p] + 273.15))); //'Bernacchi via Medlyn
      rday25 = vmax25 * 0.01; //'from Medlyn 2002
      rdaysh[p] = rday25 * pow(2, ((leaftempsh[p] - 25) / 10.0));
      rdaysh[p] = rdaysh[p] * pow((1 + exp(1.3 * (leaftempsh[p] - 55))), -1); //'high temp inhibition, collatz
      if (night == "n" && gcancsh[p] > 0) { //if// //'solve for A and ci
         if (p == 1) { //if// //'first find the mitochondrial compensation point
            ci = comp - 0.00000001; //'start at photorespiratory compensation point
            do //'loop through A-ci curve to find mitochondrial compensation point
            {
               ci = ci + 0.00000001;
               jact = ((qmax * qsh + jmax) - pow((pow((-qmax * qsh - jmax), 2) - 4 * lightcurv * qmax * qsh * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
               je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
               numerator = vmax * (ci - comp);
               denominator = ci + kc * (1 + oa / ko);
               jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
               var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
               var = var - rdaysh[p]; //'convert gross to NET A
            } while (!(var >= 0 || ci >= ca));
            //Loop Until var >= 0 Or ci >= ca //'there//'s positive Anet
            psynsh[p] = var; //'always start with zero or above
            cinsh[p] = ci; //'the dark compensation point...not predicting negative Anet
         } //End if// //'p=0 if
         if (p > 1) { //if//
            ci = cinsh[p - 1]; //'start from previous ci
            do //'loop through A-ci curve to find mitochondrial compensation point
            {
               ci = ci + 0.00000001;
               marker = gcancsh[p] * (ca - ci); //'marker is NET A from g (in umol m-2s-1, need to convert ca-ci to mole fraction), gets smaller as ci increase
               jact = ((qmax * qsh + jmax) - pow((pow((-qmax * qsh - jmax), 2) - 4 * lightcurv * qmax * qsh * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
               je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
               numerator = vmax * (ci - comp);
               denominator = ci + kc * (1 + oa / ko);
               jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
               var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
               var = var - rdaysh[p]; //'convert gross to NET A
            } while (!(var >= marker || ci >= ca));
            //Loop Until var >= marker Or ci >= ca //'when "a-ci" value just reaches the "g" value, that//'s the right ci...
            psynsh[p] = var; //'net assimilation, umol m-2 leaf area s-1
            if (psynsh[p] > psynmaxsh)
               psynmaxsh = psynsh[p];
            cinsh[p] = ci; //'store ci
         } //End if// //'p>1 if
      }
      else { //'it//'s night
         if (night == "y") {
            psynsh[p] = 0 - rdaysh[p];
            psynmaxsh = 0;
            cinsh[p] = ca; //'respiration accounted for at night
         }
         if (night == "n") {
            psynsh[p] = 0;
            psynmaxsh = 0;
            cinsh[p] = ca; //'it//'s day and p=0, g=0
         }
      } //End if// //'night if
        //'Cells(16 + p, 65) = psynsh[p]
   }

   void rhizocurves() //Generate soil E(P) global curve--only need to do this once
   {
      for (z = 1; z <= layers; z++)//z = 1 To layers
      {
         p1 = 0;
         k = 1;
         e = 0; //flow integral
         erh[z][0] = 0; //first array position is layer number
         krh[z][0] = kmaxrh[z];
         do
         {
            p2 = p1 + pinc;
            qtrapvg(p1, p2, s);
            e = e + s;
            erh[z][k] = e;
            x = p2;
            krh[z][k] = vg(x); //instantaneous k at upper limit of integration = derivative of flow integral (fundamental theorem of calculus)
            p1 = p2; //reset p1 for next increment
            k = k + 1;
            if (k == 100000)
               break; //Then Exit Do //avoid crashing for extreme vc//s
         } while (!(krh[z][k - 1] < kmin));
         //Loop Until krh(z, k - 1) < kmin
         pcritrh[z] = p2; //end of line for rhizo element z
      } //endfor// z
   } //endsub//


   void rootcurves() //generates fresh global E(P) curve for the element(erases history)
   {
      //clear arrays from earlier calls
      memset(er, 0, sizeof(er));//Erase er
      memset(kr, 0, sizeof(kr));//Erase kr
                                //do root elements
      for (z = 1; z <= layers; z++)//z = 1 To layers
      {
         kr[z][0] = ksatr[z]; //kr is instantaneous K from weibull function
                              //kminr(z) = ksatr(z)
         p1 = 0;
         k = 1;
         e = 0; //value of integral
         er[z][0] = 0; //first array position is layer number
         do
         {
            p2 = p1 + pinc;
            qtrapwbr(p1, p2, s);
            e = e + s;
            er[z][k] = e;
            x = p2;
            kr[z][k] = wbr(x); //weibull k at upper limit of integration = derivative of flow integral
            p1 = p2; //reset p1 for next increment
            k = k + 1;
            if (k == 100000)
               break;
            //If k = 100000 Then Exit Do //avoid crashing for extreme vc//s
         } while (!(kr[z][k - 1] < kmin));
         //Loop Until kr(z, k - 1) < kmin
         pcritr[z] = p2; //end of line for element z
      } //endfor// z
   } //endsub//

   void rootcurves_v() //generates fresh global E(P) curve for the element(erases history)
   {
      //clear arrays from earlier calls
      memset(er_v, 0, sizeof(er_v));//Erase er_v
      memset(kr_v, 0, sizeof(kr_v));//Erase kr_v
                                //do root elements
      for (z = 1; z <= layers; z++)//z = 1 To layers
      {
         kr_v[z][0] = ksatr[z]; //kr_v is instantaneous K from weibull function
                              //kminr(z) = ksatr(z)
         p1 = 0;
         k = 1;
         e = 0; //value of integral
         er_v[z][0] = 0; //first array position is layer number
         do
         {
            p2 = p1 + pinc;
            qtrapwbr(p1, p2, s);
            e = e + s;
            er_v[z][k] = e;
            x = p2;
            kr_v[z][k] = wbr(x); //weibull k at upper limit of integration = derivative of flow integral
            p1 = p2; //reset p1 for next increment
            k = k + 1;
            if (k == 100000)
               break;
            //If k = 100000 Then Exit Do //avoid crashing for extreme vc//s
         } while (!(kr_v[z][k - 1] < kmin));
         //Loop Until kr_v(z, k - 1) < kmin
         pcritr[z] = p2; //end of line for element z
      } //endfor// z
   } //endsub//

   void stemcurve(bool vCurve = false)
   {
      double *es_ptr = es;
      if (vCurve)
      {
         memset(es_v, 0, sizeof(es_v));
         es_ptr = es_v;
      }
      else
      {
         memset(es, 0, sizeof(es));
      }
      //memset(es, 0, sizeof(es));//Erase es //eliminate values from previous calls
      p1 = 0;
      k = 1;
      e = 0; //value of integral
      es_ptr[0] = 0;
      do
      {
         p2 = p1 + pinc;
         qtrapwbs(p1, p2, s);
         e = e + s;
         es_ptr[k] = e;
         x = p2;
         ksh = wbs(x); //weibull k
         p1 = p2; //reset p1 for next increment
         k = k + 1;
         if (k == 100000)
            break;
         //If k = 100000 Then Exit Do //avoid crashing for extreme vc//s
      } while (!(ksh < kmin));
      //Loop Until ksh < kmin
      pcrits = p2;//end of line for element z
   } //endsub//

   void leafcurve(bool vCurve = false)
   {
      double *el_ptr = el;
      if (vCurve)
      {
         memset(el_v, 0, sizeof(el_v));
         el_ptr = el_v;
      }
      else
      {
         memset(el, 0, sizeof(el));
      }
      //memset(el_ptr, 0, sizeof(el_ptr));//Erase el_ptr //eliminate values from previous calls
      p1 = 0;
      k = 1;
      e = 0; //value of integral
      el_ptr[0] = 0;
      do
      {
         p2 = p1 + pinc;
         qtrapwbl(p1, p2, s);
         e = e + s;
         el_ptr[k] = e;
         x = p2;
         ksh = wbl(x); //weibull k
         p1 = p2; //reset p1 for next increment
         k = k + 1;
         if (k == 100000)
            break;
         //If k = 100000 Then Exit Do //avoid crashing for extreme vc//s
      } while (!(ksh < kmin));
      //Loop Until ksh < kmin
      pcritl = p2; //end of line for element z
   } //endsub//

   void componentpcrits() //'gets pcrits
   {
      rhizocurves();
      rootcurves(); //'gets root element curves
      stemcurve(); //'gets stem element curve
      pcrits;
      leafcurve(); //'gets leaf element curve
      pcritl;

      rootcurves_v(); //erases history for md solution
      stemcurve(true);
      leafcurve(true);

      memset(ter, 0, sizeof(ter));
      memset(tkr, 0, sizeof(tkr));
      memset(tes, 0, sizeof(tes));
      memset(tel, 0, sizeof(tel));
   }


   void rhizoflow() //'gets flow through rhizosphere using global E(P) curve
   {
      plow = int(p1 / pinc); //'pressure index below target
      k = int(p1 / pinc);
      elow = erh[z][k]; //'e below target
      klow = krh[z][k];
      ehigh = erh[z][k + 1]; //'e above target
      khigh = krh[z][k + 1];
      plow = plow * pinc; //'convert index to pressure
      estart = (p1 - plow) / pinc * (ehigh - elow) + elow; //'linear interpolation of starting e
      klower = (p1 - plow) / pinc * (khigh - klow) + klow; //'linear interpolation of K(P)at lower limit of integration
      plow = int(p2 / pinc); //'pressure index below target
      k = int(p2 / pinc);
      elow = erh[z][k]; //'e below target
      klow = krh[z][k];
      ehigh = erh[z][k + 1]; //'e above target
      khigh = krh[z][k + 1];
      plow = plow * pinc; //'convert index to pressure
      efinish = (p2 - plow) / pinc * (ehigh - elow) + elow; //'linear interpolation of finishing e
      kupper = (p2 - plow) / pinc * (khigh - klow) + klow; //'linear interpolation of K(P) at upper limit of flow integration
      flow = efinish - estart; //'e upstream flow
   }

   void rootflow() //'gets flow through root using global E(P) curve
   {
      plow = int(p1 / pinc); //'pressure index below target
      k = int(p1 / pinc);
      elow = er[z][k]; //'e below target
      klow = kr[z][k];
      ehigh = er[z][k + 1]; //'e above target
      khigh = kr[z][k + 1];
      plow = plow * pinc; //'convert index to pressure
      estart = (p1 - plow) / pinc * (ehigh - elow) + elow; //'linear interpolation of starting e
      klower = (p1 - plow) / pinc * (khigh - klow) + klow; //'linear interpolation of starting K(P)(lower limit of integration)
      plow = int(p2 / pinc); //'pressure index below target
      k = int(p2 / pinc);
      elow = er[z][k]; //'e below target
      klow = kr[z][k];
      ehigh = er[z][k + 1]; //'e above target
      khigh = kr[z][k + 1];
      plow = plow * pinc; //'convert index to pressure
      efinish = (p2 - plow) / pinc * (ehigh - elow) + elow; //'linear interpolation of finishing e
      kupper = (p2 - plow) / pinc * (khigh - klow) + klow; //'linear interpolation of ending K(P)(upper limit of integration)
      flow = efinish - estart; //'e upstream flow
   }

   void compositecurve() //'stores composite E(P)curve and the element conductances
   {
      elayer[0][p] = 0; //'no root mediated flow in topmost layer
      for (z = 1; z <= layers; z++)//z = 1 To layers
      {
         if (layer[z] == 0) { //if//
            prhizo[z][p] = prh[z];
            p1 = prh[z];
            p2 = pr;
            rootflow();
            elayer[z][p] = flow; //'flow through layer
            if (flow != 0)
               kroot[z][p] = std::abs(elayer[z][p] / (pr - prhizo[z][p]));
            if (flow == 0) { //if//
               if (refilling == "y") { //if// //'for refilling, starting point is always weibull
                  x = pd[z];
                  kroot[z][p] = wbr(x);
               } //End if//
               if (refilling == "n")
                  kroot[z][p] = kminroot[z];
            } //End if//
         } //End if//
         if (layer[z] == 1) { //if//
            elayer[z][p] = 0; //'no flow
            if (layerfailure[z] == "root") { //if// //'root element has failed
               kroot[z][p] = 0; //'total cavitation in root
               prhizo[z][p] = pd[z]; //'rhizosphere pressure returns to the predawn value
            } //End if//
            if (layerfailure[z] == "rhizosphere") { //if// //'rhizosphere element has failed
               x = pr;
               kroot[z][p] = wbr(x); //'root element conductance = instantaneous conductance from weibull curve at pr
               prhizo[z][p] = pcritrh[z];
            } //End if//
         } //End if//
      } //for//z
      proot[p] = pr;
      pstem[p] = ps;
      pleaf[p] = pl;
      if (e > 0) { //if//
         kleaf[p] = e / (pl - ps); //'leaf element conductance
         kstem[p] = e / (ps - pr - pgrav); //'stem element conductance subtracting extra gravity drop
         kplant[p] = e / (pl - pleaf[0] - pgrav); //'whole plant k, subtracting extra gravity drop
      }
      else {
         if (refilling == "n") { //if//
            kleaf[p] = kminleaf;
            kstem[p] = kminstem;
            kplant[p] = kminplant;
         } //End if//
         if (refilling == "y") { //if//
            x = pl;
            kleaf[p] = wbl(x);
            x = ps;
            kstem[p] = wbs(x);
            //'kplant[p]=??? ignore kplant setting for e=0...probably not important
         } //End if//
           //'dedp[p] = kminplant
      } //End if//
      eplant[p] = e; //'total flow
      if (p > 0) { //if//
         if (pleaf[p] - pleaf[p - 1] == 0) { //if//
            test = 1;
         }
         else {
            dedp[p] = einc / (pleaf[p] - pleaf[p - 1]); //'dedp=instantaneous K of system
            dedpf[p] = dedp[p] / dedp[1]; //'fractional canopy conductance
         } //End if//
      } //End if//
      pcritsystem = pl;
      ecritsystem = e;
      total = p;
   }

   void leaftempsmd()  //'gets virgin sun layer leaf temp and leaf-to-air vpd from E, Campbell and Norman book
   {
      rabs = 0.5 * (0.5 * ssun + 0.5 * sref) + emiss * (0.5 * la + 0.5 * lg); //'total absorbed radiation for sun leaves; CN 11.14
      lambda = -42.9143 * airtemp + 45064.3; //'heat of vaporization for water at air temp in J mol-1
      grad = 0.1579 + 0.0017 * airtemp + 0.00000717 * pow(airtemp, 2); //'radiative conductance (long wave) at air temp in mol m-2 s-1
      gha = 1.4 * 0.135 * pow((wind / leafwidth), 0.5); //'heat conductance in mol m-2s-1
      emd = e * (1 / laperba) * (1.0 / 3600.0) * 55.4; //'convert to E per leaf area in mol m-2s-1
      numerator = rabs - emiss * sbc * pow((airtemp + 273.2), 4) - lambda * emd / 2.0; //'divide E by 2 because energy balance is two sided.
      denominator = sha * (grad + gha);
      leaftmd = airtemp + numerator / denominator; //'leaf temp for supply function
      lavpdmd = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
      lavpdmd = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * leaftmd)) - lavpdmd + vpd; //'leaf-to-air vpd
      if (lavpdmd < 0)
         lavpdmd = 0; //'don't allow negative lavpd
                      //'eplantl(p) = eplantl(p) * 1000 //'convert to mmol m-2 s-1
   }

   void leaftempsshademd()  //'gets virgin shade layer leaf temp and leaf-to-air vpd from E, Campbell and Norman book
   {
      rabs = 0.5 * (0.5 * sshade + 0.5 * sref) + emiss * lg; //'total absorbed radiation for shaded leaves
                                                             //'lambda = -42.9143 * airtemp + 45064.3 //'heat of vaporization for water at air temp in J mol-1
                                                             //'grad = 0.1579 + 0.0017 * airtemp + 0.00000717 * airtemp ^ 2 //'radiative conductance (long wave) at air temp in mol m-2 s-1
                                                             //'gha = 1.4 * 0.135 * (wind / leafwidth) ^ 0.5 //'heat conductance in mol m-2s-1
                                                             //'eplantl(p) = eplant(p) * (1 / laperba) * (1 / 3600) * 55.4 //'convert to E per leaf area in mol m-2s-1
      numerator = rabs - emiss * sbc * pow((airtemp + 273.2), 4) - lambda * emd / 2.0; //'divide E by 2 because energy balance is two sided.
      denominator = sha * (grad + gha);
      leaftshmd = airtemp + numerator / denominator; //'leaf temp for supply function
      lavpdshmd = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
      lavpdshmd = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * leaftshmd)) - lavpdshmd + vpd; //'leaf-to-air vpd
      if (lavpdshmd < 0)
         lavpdshmd = 0; //'don't allow negative lavpd
      emd = emd * 1000; //'convert to mmol m-2 s-1
   }

   void assimilationmd() //'gets virgin assimilation for sun leaves
   {
      //'get g from D and e
      if (lavpdmd == 0)
         gcanwmd = gmax; //'maximum g if no lavpd
      else
         gcanwmd = emd / lavpdmd; //'gcanopy in mmol m-2s-1 (leaf area)

                                  //'gcanw(p) = eplantl(p) / lavpd(p) //'gcanopy in mmol m-2s-1 (leaf area)
      gcancmd = (gcanwmd / 1.6) * 1000; //'convert to CO2 conductance in umol m-2 s-1
                                        //'adjust photosynthetic inputs for Tleaf
      comp = comp25 * exp((37830 * ((leaftmd + 273.15) - 298.15)) / (298.15 * gas * (leaftmd + 273.15))); //'Bernacchi via Medlyn
      numerator = (1 + exp((svvmax * 298.2 - hdvmax) / (gas * 298.2))) * exp((havmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftmd)));
      denominator = 1 + exp((svvmax * (273.2 + leaftmd) - hdvmax) / (gas * (273.2 + leaftmd)));
      vmax = vmax25 * numerator / denominator; //'vmax corrected via Leunig 2002
      numerator = (1 + exp((svjmax * 298.2 - hdjmax) / (gas * 298.2))) * exp((hajmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftmd)));
      denominator = 1 + exp((svjmax * (273.2 + leaftmd) - hdjmax) / (gas * (273.2 + leaftmd)));
      jmax = jmax25 * numerator / denominator; //'jmax corrected via Leunig 2002
      kc = kc25 * exp((79430 * ((leaftmd + 273.15) - 298.15)) / (298.15 * gas * (leaftmd + 273.15))); //'Bernacchi via Medlyn
      ko = ko25 * exp((36380 * ((leaftmd + 273.15) - 298.15)) / (298.15 * gas * (leaftmd + 273.15)));//'Bernacchi via Medlyn
      rday25 = vmax25 * 0.01; //'from Medlyn 2002
      rdaymd = rday25 * pow(2, ((leaftmd - 25) / 10.0));
      rdaymd = rdaymd * pow((1 + exp(1.3 * (leaftmd - 55))), -1); //'high temp inhibition, collatz
      if (night == "n" && gcancmd > 0)
      { //'solve for A and ci
         if (p == 1)  //'stomata have just opened: find the mitochondrial compensation point
         {
            ci = comp - 0.00000001; //'start at photorespiratory compensation point
            do //'loop through A-ci curve to find mitochondrial compensation point
            {
               ci = ci + 0.00000001;
               jact = ((qmax * qsl + jmax) - pow((pow((-qmax * qsl - jmax), 2) - 4 * lightcurv * qmax * qsl * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
               je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
               numerator = vmax * (ci - comp);
               denominator = ci + kc * (1 + oa / ko);
               jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
               var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
               var = var - rdaymd; //'convert gross to NET A
            } while (!(var >= 0 || ci >= ca));
            //Loop Until var >= 0 Or ci >= ca //'there's positive Anet or not enough light
            psynmd[p] = var; //'always start with zero or above
            cinmd = ci; //'the dark compensation point...not predicting negative Anet
         } //'p=1 if
         if (p > 1)
         {
            ci = cinmd - 0.00000001; //'cin(p - 1) //'start from previous ci
            do //'loop through A-ci curve to find ci
            {
               ci = ci + 0.00000001;
               marker = gcancmd * (ca - ci); //'marker is NET A from g (in umol m-2s-1, need to convert ca-ci to mole fraction), gets smaller as ci increase
               jact = ((qmax * qsl + jmax) - pow((pow((-qmax * qsl - jmax), 2) - 4 * lightcurv * qmax * qsl * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
               je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
               numerator = vmax * (ci - comp);
               denominator = ci + kc * (1 + oa / ko);
               jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
               var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
               var = var - rdaymd; //'convert gross to NET A
            } while (!(var >= marker || ci >= ca));
            //Loop Until var >= marker Or ci >= ca //'when "a-ci" value just reaches the "g" value, that's the right ci...
            psynmd[p] = var; //'net assimilation, umol m-2 leaf area s-1
            if (psynmd[p] > psynmaxmd)
               psynmaxmd = psynmd[p];
            cinmd = ci; //'store ci
         } //'p>1 if
      }
      else //'it's night or stomata are closed
      {
         if (night == "y")
         {
            psynmd[p] = 0 - rdaymd;
            psynmaxmd = 0; //'respiration accounted for at night
         }
         if (night == "n")
         {
            psynmd[p] = 0;
            psynmaxmd = 0; //'it's day and p=0, g=0
         }
      } //'night if
        //'Cells(16 + p, 65) = psynsh(p)
   }


   void assimilationshademd() //'gets virgin assimilation for shade leaves for midday solution
   {
      //'get g from D and e
      if (lavpdshmd == 0)
      {
         gcanwshmd = gmax; //'set to maximum if no vpd
      }
      else
      {
         gcanwshmd = emd / lavpdshmd; //'gcanopy in mmol m-2s-1 (leaf area)
      }
      //'gcanwsh(p) = eplantl(p) / lavpdsh(p) //'gcanopy in mmol m-2s-1 (leaf area)
      gcancshmd = (gcanwshmd / 1.6) * 1000; //'convert to CO2 conductance in umol m-2 s-1
                                            //'adjust photosynthetic inputs for Tleaf
      comp = comp25 * exp((37830 * ((leaftshmd + 273.15) - 298.15)) / (298.15 * gas * (leaftshmd + 273.15))); //'Bernacchi via Medlyn
      numerator = (1 + exp((svvmax * 298.2 - hdvmax) / (gas * 298.2))) * exp((havmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftshmd)));
      denominator = 1 + exp((svvmax * (273.2 + leaftshmd) - hdvmax) / (gas * (273.2 + leaftshmd)));
      vmax = vmax25 * numerator / denominator; //'vmax corrected via Leunig 2002
      numerator = (1 + exp((svjmax * 298.2 - hdjmax) / (gas * 298.2))) * exp((hajmax / (gas * 298.2)) * (1 - 298.2 / (273.2 + leaftshmd)));
      denominator = 1 + exp((svjmax * (273.2 + leaftshmd) - hdjmax) / (gas * (273.2 + leaftshmd)));
      jmax = jmax25 * numerator / denominator; //'jmax corrected via Leunig 2002
      kc = kc25 * exp((79430 * ((leaftshmd + 273.15) - 298.15)) / (298.15 * gas * (leaftshmd + 273.15))); //'Bernacchi via Medlyn
      ko = ko25 * exp((36380 * ((leaftshmd + 273.15) - 298.15)) / (298.15 * gas * (leaftshmd + 273.15))); //'Bernacchi via Medlyn
      rday25 = vmax25 * 0.01; //'from Medlyn 2002
      rdayshmd = rday25 * pow(2, ((leaftshmd - 25) / 10.0));
      rdayshmd = rdayshmd * pow((1 + exp(1.3 * (leaftshmd - 55))), -1); //'high temp inhibition, collatz
      if (night == "n" && gcancshmd > 0)
      { //'solve for A and ci
         if (p == 1) { //'first find the mitochondrial compensation point
            ci = comp - 0.00000001; //'start at photorespiratory compensation point
            do //'loop through A-ci curve to find mitochondrial compensation point
            {
               ci = ci + 0.00000001;
               jact = ((qmax * qsh + jmax) - pow((pow((-qmax * qsh - jmax), 2) - 4 * lightcurv * qmax * qsh * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
               je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
               numerator = vmax * (ci - comp);
               denominator = ci + kc * (1 + oa / ko);
               jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
               var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
               var = var - rdayshmd; //'convert gross to NET A
            } while (!(var >= 0 || ci >= ca));
            //Loop Until var >= 0 Or ci >= ca //'there's positive Anet
            psynshmd[p] = var; //'always start with zero or above
            cinshmd = ci; //'the dark compensation point...not predicting negative Anet
         } //'p=1 if
         if (p > 1)
         {
            ci = cinshmd - 0.00000001; //'start from previous ci backed off a bit
            do //'loop through A-ci curve to find mitochondrial compensation point
            {
               ci = ci + 0.00000001;
               marker = gcancshmd * (ca - ci); //'marker is NET A from g (in umol m-2s-1, need to convert ca-ci to mole fraction), gets smaller as ci increase
               jact = ((qmax * qsh + jmax) - pow((pow((-qmax * qsh - jmax), 2) - 4 * lightcurv * qmax * qsh * jmax), 0.5)) / (2 * lightcurv); //'0.9 is curvature of light response curve...this from Medlyn 2002
               je = (jact * (ci - comp)) / (4 * (ci + 2 * comp)); //'light limited PS
               numerator = vmax * (ci - comp);
               denominator = ci + kc * (1 + oa / ko);
               jc = numerator / denominator; //'rubisco limited A, umol m-2 s-1
               var = (je + jc - pow((pow((je + jc), 2) - 4 * thetac * je * jc), 0.5)) / (2 * thetac); //'gross photosynthetic rate, gets larger with ci
               var = var - rdayshmd; //'convert gross to NET A
            } while (!(var >= marker || ci >= ca));
            //Loop Until var >= marker Or ci >= ca //'when "a-ci" value just reaches the "g" value, that's the right ci...
            psynshmd[p] = var; //'net assimilation, umol m-2 leaf area s-1
            if (psynshmd[p] > psynmaxshmd)
               psynmaxshmd = psynshmd[p];
            cinshmd = ci; //'store ci
         } //'p>1 if
      }
      else //'it's night
      {
         if (night == "y")
         {
            psynshmd[p] = 0 - rdayshmd;
            psynmaxshmd = 0; //'respiration accounted for at night
         }
         if (night == "n")
         {
            psynshmd[p] = 0;
            psynmaxshmd = 0; //'it's day and p=0, g=0
         }
      } //'night if
        //'Cells(16 + p, 65) = psynsh(p)
   }



   void ludcmp() //'does LU decomposition on the jacobian prior to solution by lubksb
   {
      d = 1;
      for (i = 1; i <= unknowns; i++)//i = 1 To unknowns
      {
         aamax = 0;
         for (j = 1; j <= unknowns; j++)//j = 1 To unknowns
         {
            if (std::abs(jmatrix[i][j]) > aamax)
               aamax = std::abs(jmatrix[i][j]);
         }
         if (aamax == 0)
            return;
         vv[i] = 1 / aamax;
      }
      for (j = 1; j <= unknowns; j++)//j = 1 To unknowns
      {
         for (i = 1; i < j; i++) //i = 1 To j - 1
         {
            sum = jmatrix[i][j];
            for (k = 1; k < i; k++) //k = 1 To i - 1
            {
               sum = sum - jmatrix[i][k] * jmatrix[k][j];
            }
            jmatrix[i][j] = sum;
         }
         aamax = 0;
         for (i = j; i <= unknowns; i++)//i = j To unknowns
         {
            sum = jmatrix[i][j];
            for (k = 1; k < j; k++) //k = 1 To j - 1
            {
               sum = sum - jmatrix[i][k] * jmatrix[k][j];
            }
            jmatrix[i][j] = sum;
            dum = vv[i] * std::abs(sum);
            if (dum > aamax)
            {
               imax = i;
               aamax = dum;
            }
         }
         if (j != imax) // j <> imax
         {
            for (k = 1; k <= unknowns; k++) //k = 1 To unknowns
            {
               dum = jmatrix[imax][k];
               jmatrix[imax][k] = jmatrix[j][k];
               jmatrix[j][k] = dum;
            }
            d = -d;
            vv[imax] = vv[j];
         }
         indx[j] = imax;
         if (jmatrix[j][j] == 0)
            jmatrix[j][j] = 1E-25;
         if (j != unknowns)
         {
            dum = 1 / jmatrix[j][j];
            for (i = j + 1; i <= unknowns; i++) //i = j + 1 To unknowns
            {
               jmatrix[i][j] = jmatrix[i][j] * dum;
            }
         }
      }
      memset(vv, 0, sizeof(vv));
      //Erase vv
   }

   void lubksb() //'solves the decomposed jacobian for delta p's
   {
      ii = 0;
      for (i = 1; i <= unknowns; i++)//i = 1 To unknowns
      {
         ll = int(indx[i]); //'indx array comes from ludcmp
         sum = func[ll]; //'the func array input is the right-hand vector
         func[ll] = func[i];
         if (ii != 0)
         {
            for (j = ii; j < i; j++)//j = ii To i - 1
            {
               sum = sum - jmatrix[i][j] * func[j];
            }
         }
         else
         {
            if (sum != 0)
               ii = i;
         }
         func[i] = sum;
      }
      for (i = unknowns; i >= 1; i--) //i = unknowns To 1 Step -1
      {
         sum = func[i];
         for (j = i + 1; j <= unknowns; j++) //j = i + 1 To unknowns
         {
            sum = sum - jmatrix[i][j] * func[j];
         }
         func[i] = sum / jmatrix[i][i];
      }
   }

   void newtonrhapson() //returns rhizosphere pressures and root pressure, pr, as function of pd's and e
   {

      //'prinitial = pr //record the original guess
      weird = 0; //tracks pr estimate
      check = 0; //restart loop counter
      do //loop to reset guesses
      {
         //'restore layer functioning
         for (z = 1; z <= layers; z++)
         {
            if (kminroot[z] != 0)
            {
               layer[z] = 0; //make sure to start with all layers functioning
               layerfailure[z] = "no failure"; //reset
            }
         }
         failspot = "no failure";

         if (weird == 1) //reset guesses
         {

            double rFloat = (double)rand() / (double)RAND_MAX;
            k = int((layers - 1 + 1) * rFloat + 1);
            pr = pd[k]; //random choice of pd

            if (false) // can enable this is running into erroneous solutions -- but allowing these to vary instead of resetting results in more frequent solutions in my experience
            { // alternatively could randomize them properly
               for (k = 1; k <= layers; k++) //reset prhz(k)
               {
                  prh[k] = pd[k];
               }
            }

            weird = 0; //reset cutoff
         } //end reset guesses loop
         check = check + 1; //number of restarts
         ticks = 0; //convergence counter
         do //loop to seek convergence
         {
            ticks = ticks + 1;

            if (ticks > 1000)
            {
               weird = 1;
               std::cout << "NR ticks exceeded 1000 -- setting weird=1 for retry. Pinc too high? dd = " << dd << std::endl;
            }
            //'get top row of matrix and right-hand func vector
            //'zero out the jacobian first
            for (k = 1; k <= unknowns; k++)
            {
               for (j = 1; j <= unknowns; j++)
                  jmatrix[k][j] = 0;
            }

            //'fill up four arrays:
            //'func(i) is zero flow function...the right-hand flow vector
            //'dfrhdprh(i) is partial derivative of frh(i) for prh(i)(rhizo pressure)...the diagonal of the jacobian
            //'dfrhdpr(i) is partial derivative of frh(i)for pr (root pressure)...the last column of the jacobian
            //'dfrdprh(i) is partial derivative of fr for prh(i)...the last row of the jacobian
            frt = 0; //this is the last row of right-hand flow vector
            dfrdpr = 0; //this is lower right-hand partial for jacobian
            for (z = 1; z <= layers; z++)
            {
               if (pd[z] >= pcritrh[z] || prh[z] >= pcritrh[z])
               {
                  layer[z] = 1; //layer's just gone out of function
                  layerfailure[z] = "rhizosphere";
                  //weird = 1; //could be result of non-convergence
               }
               if (layer[z] == 0)  //it's functional
               {
                  p1 = pd[z]; //p1 is not a guess
                  p2 = prh[z]; //prh[z] IS a guess and it is initially set prior to the RN routine
                  rhizoflow(); //gets flows through rhizosphere element from p2-p1 and E(p)curve; gets K's at p1 and p2 as well
                  func[z] = flow;
                  dfrhdprh[z] = kupper;
               }
               if (prh[z] >= pcritr[z] || pr >= pcritr[z])
               {
                  layer[z] = 1; //layer's just gone out of function
                  layerfailure[z] = "root";
                  //weird = 1; //could be result of non-convergence
               }
               if (layer[z] == 0) //it's functional
               {
                  p1 = prh[z]; //now re-set p1 to prh[z]...the guess
                  p2 = pr; //guess must come from before NR routine?
                  rootflow(); //gets flows through root element, and K's at upper and lower bounds
                  func[z] = func[z] - flow;
                  dfrhdprh[z] = dfrhdprh[z] + klower;
                  dfrhdpr[z] = -kupper;
                  dfrdprh[z] = -klower;
               }
               if (layer[z] == 1)  //layer's out of function...zero out values
               {
                  func[z] = 0;
                  dfrhdprh[z] = 1E-25;
                  dfrdprh[z] = 1E-25;
                  dfrhdpr[z] = 1E-25;
                  kupper = 1E-25;
                  flow = 0;
               }
               dfrdpr = dfrdpr + kupper;
               frt = frt + flow;
            }
            frt = frt - e;
            //'now load jacobian
            for (k = 1; k <= layers; k++)
            {
               jmatrix[k][unknowns] = dfrhdpr[k]; //last column with dFrh/dPr partials
            }
            for (k = 1; k <= layers; k++)
            {
               jmatrix[unknowns][k] = dfrdprh[k]; //last row with dFr/dPrh partials
            }
            for (k = 1; k <= layers; k++)
            {
               jmatrix[k][k] = dfrhdprh[k]; //diagonal of dFrh/dPrh partials
            }
            jmatrix[unknowns][unknowns] = dfrdpr; //lower right corner with dFr/dPr partial
            func[unknowns] = frt; //last position in right-hand flow vector

                                  //'ok, jacobian and righthand vector are loaded
                                  //'test for total failure
            sum = 0;
            for (k = 1; k <= layers; k++)
            {
               sum = sum + layer[k];
            }
            if (sum == layers)  //total failure
            {
               failspot = "belowground";
               weird = 1; //trigger a restart
            }
            //'test for flow conservation (steady-state)
            threshold = 0;
            for (k = 1; k <= unknowns; k++) //k = 1 To unknowns
            {
               threshold = threshold + std::abs(func[k]);
            }
            if (ticks == 1)
               initialthreshold = threshold;
            //'remember to replace "n" with "unknowns" in ludcmp and lubksb
            ludcmp(); //numerical recipe for doing LU decomposition of jacobian prior to solving
            lubksb(); //solves the decomposed jacobian for delta p's
                      //'print out solution vector of pressures
                      //'revise unknown pressures
            for (k = 1; k <= layers; k++)//k = 1 To layers
            {
               prh[k] = prh[k] - func[k]; //NOTE lubksb replaces original right-side func()vector with the solution vector
            }
            pr = pr - func[unknowns];
            //'check for jumping lower bound
            for (k = 1; k <= layers; k++)//k = 1 To layers
            {
               if (prh[k] < 0)
                  prh[k] = 0;
            }
            if (pr < 0)
               pr = 0;
            //'if pr > pcritr Then
            //'pr = prinitial
            //'weird = 1 //trigger a re-start
            //'}
            if (ticks > 1)  //check for non convergence
            {
               if (threshold > initialthreshold)
                  weird = 1;//pr is spiraling, restart NR with new guesses
               if (pr >= pcrits)
                  weird = 1; //
            }
         } while (!(threshold < 0.02 || weird == 1));
         //Loop Until threshold < 0.01 Or weird = 1 //weird = 1 restarts NR


      } while (!((threshold < 0.02 && weird == 0) || check > 500));

      if (check > 500)
      {
         // disable this output if it's causing too much spam
         //std::cout << "NR Failure " << threshold << " check = " << check << " dd = " << dd << " watercontent = " << waterold * 1000.0 << " weird = " << weird << std::endl;
         // keep track of the frequency of NR failures
         gs_ar_nrFailConverge[gs_yearIndex]++; // non-convergent failure
         gs_ar_nrFailConverge_Water[gs_yearIndex] += waterold * 1000.0;
         if (waterold * 1000.0 > gs_ar_nrFailConverge_WaterMax[gs_yearIndex])
            gs_ar_nrFailConverge_WaterMax[gs_yearIndex] = waterold * 1000.0;

      }

      //Loop Until threshold < 0.01 And weird = 0 Or check > 500 //give up after 2000 restarts
      //'if check >= 500 Then Stop

      //final step -- recheck the layers
      for (z = 1; z <= layers; z++)
      {
         if (kminroot[z] != 0)
         {
            layer[z] = 0; //make sure to start with all layers functioning
            layerfailure[z] = "no failure"; //reset
         }

         if (pd[z] >= pcritrh[z] || prh[z] >= pcritrh[z])
         {
            layer[z] = 1; //layer's just gone out of function
            layerfailure[z] = "rhizosphere";
         }
         if (prh[z] >= pcritr[z] || pr >= pcritr[z])
         {
            layer[z] = 1; //layer's just gone out of function
            layerfailure[z] = "root";
         }
      }
   }

   void stem() //gets stem pressure and conductance from pr and e.
   {
      //'start with stem
      p1 = pr + pgrav; //add gravity drop before integration
      plow = int(p1 / pinc); //pressure index below target
      k = int(p1 / pinc);
      elow = es[k]; //e below target
      ehigh = es[k + 1]; //e above target
      plow = plow * pinc; //convert index to pressure
      estart = (p1 - plow) / pinc * (ehigh - elow) + elow; //linear interpolation of starting e
      efinish = estart + e;
      j = k;
      do //find efinish
      {
         j = j + 1;
         if (es[j] == 0)
         {
            test = 1;
            failspot = "stem";
            return;
         }
      } while (!(es[j] > efinish));
      //Loop Until es(j) > efinish
      ehigh = es[j];
      elow = es[j - 1];
      p2 = ((efinish - elow) / (ehigh - elow)) * pinc + pinc * (j - 1); //pstem
      ps = p2; //ps is downstream stem pressure
      if (ps >= pcrits)
      {
         test = 1;
         failspot = "stem";
         return;
      }
   }

   void leaf() //gets leaf pressure from stem pressure and e
   {
      p1 = ps;
      plow = int(p1 / pinc); //pressure index below target
      k = int(p1 / pinc);
      elow = el[k]; //e below target
      ehigh = el[k + 1]; //e above target
      plow = plow * pinc; //convert index to pressure
      estart = (p1 - plow) / pinc * (ehigh - elow) + elow; //linear interpolation of starting e
      efinish = estart + e;
      j = k;
      do //find efinish
      {
         j = j + 1;
         if (el[j] == 0) //= Empty 
         {
            test = 1;
            failspot = "leaf";
            return;
         }
      } while (!(el[j] > efinish));
      //Loop Until el(j) > efinish
      ehigh = el[j];
      elow = el[j - 1];
      p2 = ((efinish - elow) / (ehigh - elow)) * pinc + pinc * (j - 1); //pleaf
      pl = p2; //pl is leaf pressure
      if (pl >= pcritl)
      {
         test = 1;
         failspot = "leaf";
      }
   }

   void storehistory() //'stores historical element e(P) curves
   {
      /*memset(ter, 0, sizeof(ter));
      memset(tkr, 0, sizeof(tkr));
      memset(tes, 0, sizeof(tes));
      memset(tel, 0, sizeof(tel));*/

      /*memcpy(ter, er, sizeof(ter));
      memcpy(tkr, kr, sizeof(tkr));
      memcpy(tes, es, sizeof(tes));
      memcpy(tel, el, sizeof(tel));

      memcpy(er, er_v, sizeof(er));
      memcpy(kr, kr_v, sizeof(kr));
      memcpy(es, es_v, sizeof(es));
      memcpy(el, el_v, sizeof(el));*/

      for (z = 1; z <= layers; z++)
      {
         k = 0;
         do
         {
            ter[z][k] = er[z][k];
            tkr[z][k] = kr[z][k];

            er[z][k] = er_v[z][k];
            kr[z][k] = kr_v[z][k];
            k = k + 1;
            if (k == 100000)
               break;
         } while (er[z][k] != 0 || er_v[z][k] != 0 || ter[z][k] != 0);
         //Loop Until er(z, k) = Empty
      }
      k = 0;
      do
      {
         tes[k] = es[k];

         es[k] = es_v[k];
         k = k + 1;
         if (k == 100000)
            break;
      } while (es[k] != 0 || es_v[k] != 0 || tes[k] != 0);
      //Loop Until es(k) = Empty
      k = 0;
      do
      {
         tel[k] = el[k];

         el[k] = el_v[k];
         k = k + 1;
         if (k == 100000)
            break;
      } while (el[k] != 0 || el_v[k] != 0 || tel[k] != 0);
      //Loop Until el(k) = Empty
      for (z = 1; z <= layers; z++)
      {
         tlayer[z] = layer[z];
         tlayerfailure[z] = layerfailure[z];
      }
      for (z = 1; z <= layers; z++) // 'get all layers functioning to start with
      {
         if (layerfailure[z] == "ok")
            layer[z] = 0;
         else if (layerfailure[z] == "rhizosphere")
            layer[z] = 0;
         else if (layerfailure[z] == "root" && kminroot[z] == 0)
            layer[z] = 1; // 'take out root layer if failed at start of composite curve
         else
            layer[z] = 0; //'still functional
      }
   }

   void gethistory() //'restores historical element e(P) curves
   {
      for (z = 1; z <= layers; z++)
      {
         k = 0;
         do
         {
            er[z][k] = ter[z][k];
            kr[z][k] = tkr[z][k];
            k = k + 1;
            if (k == 100000)
               break;
         } while (ter[z][k] != 0 || er[z][k] != 0);
         //Loop Until ter[z][k] = Empty
      }
      k = 0;
      do
      {
         es[k] = tes[k];
         k = k + 1;
         if (k == 100000)
            break;
      } while (tes[k] != 0 || es[k] != 0);
      //Loop Until tes[k] = Empty
      k = 0;
      do
      {
         el[k] = tel[k];
         k = k + 1;
         if (k == 100000)
            break;
      } while (tel[k] != 0 || el[k] != 0);
      //Loop Until tel[k] = Empty
      //'re-set failure status (this at the critical point)
      for (z = 1; z <= layers; z++)
      {
         layer[z] = tlayer[z];
         layerfailure[z] = tlayerfailure[z];
      }
   }

   void canopypressure()
   {
      //computes carbon-based middays; pressures and cost curve from virgin supply function, gas exchange from historical values
      //store history, get MD and cost function from virgin curve
      if (ecritsystem == 0)
      {
         k = 0;
         goto tenMarker; //don't bother to find midday
      }
      storehistory(); //stores xylem element curves and failure status, resets layers to full functioning
      //rootcurves(); //erases history for md solution
      //stemcurve();
      //leafcurve();
      sum = 0;
      t = 0;
      for (k = 1; k <= layers; k++) //assign source pressures, set layer participation
      {
         if (layer[k] == 0)
         {
            prh[k] = pd[k]; //initial guess of unknown rhizosphere pressures
            sum = sum + pd[k];
         }
         else
         {
            prh[k] = pcritr[k];
            t = t + 1;
         }
      }
      if (t < layers)
      {
         pr = sum / (layers - t); //set unknown proot to average pd
      }
      else
      {
         failure = 1; //system is critical
         return;
      }
      test = 0; //=1 if stem or leaf fails
                //virgin gain function params
      psynmaxmd = 0;
      psynmaxshmd = 0;
      //now loop through virgin risk curve
      e = -einc;
      p = -1;
      dedplmin = ksatp; //insures the kloss function is monotonic
      do
      {
         e = e + einc;
         p = p + 1;
         newtonrhapson(); //pd's already assigned...this solves for p's and e's in fingers of chicken foot

         if (check >= 500)
         {
            //p = p - 1
            //std::cout << "NR failed on virgin curve at e = " << e << " on timestep dd = " << dd << std::endl;
            break; //gone as far as can go
         }
         stem(); //gets stem and leaf pressures
         leaf();
         pleafv[p] = pl; //pleaf from virgin curve
         if (p == 0)
         {
            predawn = pl; //set the predawn...pl returned by "leaf" routine
            plold = pl;
         }
         if (p > 0)
         {
            if ((pl - plold) == 0)
               break; //gone to failure
            if (p == 1)
               dedplzero = einc / (pl - plold); //note: pl is returned by "leaf" routine
            dedpl = einc / (pl - plold);
            if (dedpl < dedplmin)
               dedplmin = dedpl; //insure that kloss only goes up
            klossv[p] = dedplzero - dedpl; //non-normalized kloss from virgin curve
                                           //dedpf(p) = dedpl / dedplzero //fractional k/kmax canopy from virgin curve (units don't matter!)
         }
         if (pl >= pcritsystem)
            break; //gone to failure
                   //now get virgin A curve

         leaftempsmd(); //gets virgin sun layer leaf temperature from energy balance
         leaftempsshademd(); //gets virgin shade layer leaf temperature
         assimilationmd(); //gets virgin sun layer photosynthesis
         assimilationshademd(); //gets virgin shade layer photosynthesis
                                //by now we have assigned psynmd(p) and psynshmd(p) and reset psynmaxes

         plold = pl;
         //e = e + einc
         //p = p + 1
      } while (!(test == 1 || p > 99900)); //loop to failure or out of "p"s
                                           //Loop Until test = 1 Or p > 99900 'loop to failure or out of "p"s
                                           //If check >= 2000 Then Exit Sub

      if (p <= 2)
      {
         k = 0; //too close to critical, shut down system
         goto tenMarker;
      }

      klossv[0] = 0;
      maxkloss = dedplzero - dedpl; //maximum kloss...may not be kmax
      totalv = p - 1;
      klossv[totalv] = maxkloss;
      //now normalize klossv for virgin pleafv
      for (p = 0; p <= totalv; p++)
      {
         klossv[p] = klossv[p] / maxkloss;
      }

      //now, find the middday from virgin risk and gain
      //First do for sun layer
      p = -1; //p is still index for virgin curve
      dpmax = 0;
      dpamax = -100;
      amaxmax = 0; //insures the gain function is monotonic

      dpamin = 0.0; // keep track of low values to avoid extreme negative profit curves producing a result
      do
      {
         p = p + 1;
         amaxfrac[p] = psynmd[p] / psynmaxmd; //this is the normalized revenue function from the virgin curve
         if (amaxfrac[p] > amaxmax)
         {
            amaxmax = amaxfrac[p];
         }
         else
         {
            amaxfrac[p] = amaxmax;
         } //this insures that amaxfrac monotonically increases
         if (amaxfrac[p] < 0)
            amaxfrac[p] = 0; //no negative gains
         if (klossv[p] < 0)
            klossv[p] = 0; //no negative risks
         dpa[p] = amaxfrac[p] - klossv[p]; //profit, with revenue from virgin curve and cost from virgen one
         if (p < runmean - 1)
            rmean = 0;
         if (p >= runmean - 1)  //get running mean
         {
            sum = 0;
            for (i = p - runmean + 1; i <= p; i++)
            {
               sum = sum + dpa[i];
            }
            rmean = sum / runmean; //the running mean
         }
         if (rmean < dpamin)
         {
            dpamin = rmean;
         }
         if (rmean < 0)
            rmean = 0; //avoid negative rmean
         if (rmean > dpamax)
         {
            dpamax = rmean;
            md = pleafv[p]; //midday pressure for sun layer from virgin curves
         }
         //print out gain and cost

      } while (!(einc * p >= gmax * lavpd[p] || total == 0 || (rmean < dpamax / cutoff && p > runmean && p > 15) || klossv[p] > 0.9 || p >= totalv));

      //std::cout << "DPA MIN = " << dpamin << " DPA MAX = " << dpamax << std::endl;
      if (dpamin < 0.0 && dpamax > 0.0 && std::abs(dpamin) > dpamax)
      {
         // the profit went more negative than positive, so reset mid-day to predawn
         md = pleafv[0];
      }

      // [HNT] debug
      if (!(rmean < dpamax / cutoff))
      {
         //std::cout << "Terminated sun layer opt without finding peak! At timestep dd = " << dd << std::endl;
         if (einc * p >= gmax * lavpd[p])
         {
            ;// std::cout << "Terminated sun layer opt without finding peak: end case 1 einc * p >= gmax * lavpd[p]" << std::endl;
         }
         if (total == 0)
         {
            std::cout << "Terminated sun layer opt without finding peak: end case 2 total == 0" << std::endl;
         }
         if (p >= totalv)
         {
            std::cout << "Terminated sun layer opt without finding peak: end case 3 exceeded end of virgin curve" << std::endl;
         }
         if (klossv[p] > 0.9)
         {
            std::cout << "Terminated sun layer opt without finding peak: end case 4 klossv[p] > 0.9" << std::endl;
         }
      }
      // [/HNT]

      //while (!(einc * p >= gmax * lavpd[p] || total == 0 || rmean < dpamax / cutoff && p > runmean || klossv[p] > 0.9 || p >= totalv));
      //Loop Until einc * p >= gmax * lavpd(p) Or total = 0 Or rmean < dpamax / cutoff And p > runmean Or klossv(p) > 0.9 Or p >= totalv //loop until g maxed out or to failure...note e and g in kg hr-1

      dpasun = dpamax; // unused?

                       //now do for shade layer
      if (psynmaxshmd == 0) //shade layer's below light compensation point, don't open
      {
         mdsh = pleafv[0]; //set midday = predawn
      }
      else
      {
         p = -1; //p is still index for virgin curve
         dpmax = 0;
         dpamax = -100;
         amaxmax = 0;

         dpamin = 0.0;
         do
         {
            p = p + 1;
            amaxfracsh[p] = psynshmd[p] / psynmaxshmd; //this is the normalized shade revenue function from the virgin curve
            if (amaxfracsh[p] > amaxmax)
            {
               amaxmax = amaxfracsh[p];
            }
            else
            {
               amaxfracsh[p] = amaxmax;
            }  //this insures that amaxfrac monotonically increases
               //now loop to find kloss from virgin curve that matches historical pleaf
            if (amaxfracsh[p] < 0)
               amaxfracsh[p] = 0; //no negative gains
            if (klossv[p] < 0)
               klossv[p] = 0; //no negative risks
            dpa[p] = amaxfracsh[p] - klossv[p]; //profit, with revenue from historical curve and cost from virgen one
            if (p < runmean - 1)
               rmean = 0;
            if (p >= runmean - 1)  //get running mean
            {
               sum = 0;
               for (i = p - runmean + 1; i <= p; i++)
               {
                  sum = sum + dpa[i];
               }
               rmean = sum / runmean; //the running mean
            }
            if (rmean < dpamin)
            {
               dpamin = rmean;
            }
            if (rmean < 0)
               rmean = 0; //avoid negative dpa
            if (rmean > dpamax)
            {
               dpamax = rmean;
               mdsh = pleafv[p]; //midday pressure for shade layer from virgin curve
            }

         } while (!(einc * p >= gmax * lavpdsh[p] || total == 0 || (rmean < dpamax / cutoff && p > runmean && p > 15) || klossv[p] > 0.9 || p >= totalv));
         //Loop Until einc * p >= gmax * lavpdsh[p] Or total = 0 Or rmean < dpamax / cutoff And p > runmean Or klossv[p] > 0.9 Or p >= totalv //loop until g maxed out or to failure...note e and g in kg hr-1
         if (dpamin < 0.0 && dpamax > 0.0 && std::abs(dpamin) > dpamax)
         {
            // the profit went more negative than positive, so reset mid-day to predawn
            mdsh = pleafv[0];
         }
      } //psynmaxsh if

      k = -1;
      //Range("c17:f10000").ClearContents
      do
      {
         k = k + 1;
      } while (!(pleaf[k] >= md || pleaf[k] == 0));
      //Loop Until pleaf(k) >= md Or pleaf(k) = Empty //pleaf from historical curve must match pleaf from virgin curve
      transpiration = eplantl[k]; //all gas exchange values are from most recent historical values
      psynact = psyn[k];
      gcmd = gcanw[k]; //g for water in mmol
      lavpdmd = lavpd[k] * patm;
      cinc = cin[k];
      //If k > 1 Then deda = (eplantl(k) - eplantl(k - 1)) / (psyn(k) - psyn(k - 1))
      halt = k; //halt is index of midday datum
                //now do shade layer
      k = -1;
      do
      {
         k = k + 1;
      } while (!(pleaf[k] >= mdsh || pleaf[k] == 0));
      //Loop Until pleaf(k) >= mdsh Or pleaf(k) = Empty //pleaf for historical curve must match pleaf from virgin curve
      transpirationsh = eplantl[k]; //all gas exchange values are from most recent historical values
      psynactsh = psynsh[k];
      gcmdsh = gcanwsh[k]; //g for water in mmol
      lavpdmdsh = lavpdsh[k] * patm;
      cincsh = cinsh[k];
      //If k > 1 Then deda = (eplantl(k) - eplantl(k - 1)) / (psyn(k) - psyn(k - 1))
      haltsh = k; //halt is index of midday datum
      if (ecritsystem == 0)
      {
      tenMarker:     //no midday
         k = 0;
         transpiration = eplantl[k]; //all gas exchange values are from most recent historical values
         psynact = psyn[k];
         gcmd = gcanw[k]; //g for water in mmol
         lavpdmd = lavpd[k] * patm;
         cinc = cin[k];
         halt = k;
         transpirationsh = eplantl[k]; //all gas exchange values are from most recent historical values
         psynactsh = psynsh[k];
         gcmdsh = gcanwsh[k]; //g for water in mmol
         lavpdmdsh = lavpdsh[k] * patm;
         cincsh = cinsh[k];
         haltsh = k; //halt is index of midday datum
      }
      gethistory(); //reinstates historical element curves and failure status prior to updating
      if (refilling == "y") //need to record midday kmins, uses sunlit pressures
      {
         for (z = 1; z <= layers; z++)
         {
            kminroot[z] = kroot[z][halt];
         }
         kminstem = kstem[halt];
         kminleaf = kleaf[halt];
      }
   }

   long modelTimestepIter(long &VBA_dd)
   {
      dd = VBA_dd;

      if (dd == 1 || isNewYear)
      {
         failure = 0;//'=1 for system failure at midday...terminates run
         failspot = "no failure";
         componentpcrits();//'gets pcrits for each component
         failspot = "no failure";

         for (k = 1; k <= layers; k++)// k = 1 To layers ;//'exclude the top layer
         {
            kminroot[k] = ksatr[k];
         }

         kminstem = ksats;
         kminleaf = ksatl;
         kminplant = ksatp;

         gwflow = 0; //'inflow to bottom of root zone
         drainage = 0; //'drainage from bottom of root zone

         gs_prevDay = 0;
         gs_inGrowSeason = false;
         gs_doneFirstDay = false; // prevents PLC misses on first year
                                  //[/HNT]
      }

      //dd++;
      long testCount = 0;
      yearVal = std::lround(dSheet.Cells(rowD + dd, colD + dColYear));
      if (yearVal != year_cur)
      {
         // all of these year values get zeroed out during initModelVars, so we can assume they will be zero at the start of a new iteration and test against 0 meaning "not set"
         if (yearVal > year_cur && yearVal > year_start)
         {
            // if the start year is zero, this is the first year we've processed
            if (year_start == 0) // model runs starting in the year zero are not supported, I guess?
               year_start = yearVal;
            year_cur = yearVal;

            gs_yearIndex = year_cur - year_start;
            if (gs_yearIndex >= 100)
               gs_yearIndex = 99; //safety check

            gs_ar_years[gs_yearIndex] = year_cur; // correct the year listing in the "growing seasons" array 
                                                  // TODO make that input data able to handle the new system?? Otherwise just eliminate
                                                  // on VBA side it might be sufficient to pull the start year from the first line of data
                                                  // though this still assumes that our data is in linear time order ...

            if (gs_yearIndex > 0)
            {
               // if we've set the year and it was anythign other than the first timestep of this iteration (in which case we set the year index from zero TO zero)
               // then return to VBA and let it handle the model reset
               // it will re-run this timestep afterwards
               gs_prevDay = 0;
               gs_inGrowSeason = false;
               gs_doneFirstDay = false;

               isNewYear = true; // this can be used to override the dd==/>1 cases -- it will be set to false upon successful completion of a timestep
               return gs_yearIndex;
            }
            // if we detected a year change, but the yearIndex comes out to zero, then this must be the first timestep
            // of the first year in the data -- so we can continue without returning and re-setting the model
         }
         else // going back in time means something went wrong -- may be able to handle out of order years later though
         {
            return 0; // failure
         }
      }

      jd = dSheet.Cells(rowD + dd, colD + dColDay); //'julian day

                                                    //Debug.Print "DOING A LOOP-1 " & dd

      if (dd > 1 && !isNewYear) { //if// //'set timestep
         if (tod < dSheet.Cells(rowD + dd, colD + dColTime)) { //if// //'same old day getting older
            timestep = dSheet.Cells(rowD + dd, colD + dColTime) - tod;
         }
         else {
            timestep = (24 - tod) + dSheet.Cells(rowD + dd, colD + dColTime); //'a new day has started
                                                                              //[HNT] multi-year support
                                                                              // following method of new year detection has been replaced with the method above
                                                                              // now that we are including the year as a data input
            gs_prevDay = jd;
            gs_inGrowSeason = true; // isInGrowSeasonSimple(); //it's a new day, so let's see if this is in the growing season or not
                                    //[/HNT]
         } //End if// //'tod if
      } //End if// //'dd>1 if
        //[HNT] multi-year support
      else
      {
         gs_inGrowSeason = true; // isInGrowSeasonSimple(); //it's the first data point, test if we're starting in growing season
      }

      gs_inGrowSeason = isInGrowSeasonSimple(); // just always call this!
                                                //[/HNT]

      tod = dSheet.Cells(rowD + dd, colD + dColTime); //'time of day, standard local time in hour fraction
      obssolar = dSheet.Cells(rowD + dd, colD + dColSolar); //'observed total solar, horizontal, wm-2
      vpd = dSheet.Cells(rowD + dd, colD + dColD); //'midday vpd in kPa
      vpd = vpd / patm; //'vpd in mole fraction
      airtemp = dSheet.Cells(rowD + dd, colD + dColTAir); //'in C
      maxvpd = (101.3 / patm) * (-0.0043 + 0.01 * exp(0.0511 * airtemp)); //'saturated mole fraction of vapor in air
      wind = dSheet.Cells(rowD + dd, colD + dColWind); //'wind speed
      if (wind < minwind) { //if//
         dSheet.Cells(rowD + dd, colD + dColWind) = minwind; //'set to minimum wind
         wind = minwind;
      } //End if//
      us = wind * 0.1; //'understory windspeed in m s-1
      soiltemp = dSheet.Cells(rowD + dd, colD + dColTSoil); //'surface temp of soil
      if (vpd > maxvpd) { //if//
         vpd = maxvpd;
         dSheet.Cells(rowD + dd, colD + dColD) = maxvpd * patm; //'print out maximum vpd
      } //End if//

      getsoilwetness(); //'after initializing, start updating water contents of soil layers
      solarcalc(); //'get radiation for timestep
      if (qsl > lightcomp /*[HNT] multiyear*/ && gs_inGrowSeason /*[/HNT]*/) { //if//
         night = "n"; //'it//'s light enough for sun layer to do business
      }
      else {
         night = "y"; //'too dark
      } //End if// //'end night if

      gwflow = 0; //'re-set inflow to bottom of root zone
      drainage = 0; //'re-set drainage from bottom of root zone
                    //'Call getpredawns //'update soil pressure of each layer
                    //'if failure = 1 { //if// Exit do
      chalk = 0; //'einc counter

   twentyMarker:

      getpredawns(); //'passed initializing...update soil pressure of each layer

      if (failure == 1)
         return -1;//break;//Exit do

      test = 0; //'=1 if stem or leaf fails
      p = -1; //'E(Pleaf) point counter
              //'set initial guesses of the three unknowns
      e = -einc; //'total e: einc and e are still in kg hr-1 m-2 basal area, as are conductances
                 //'Range("c18:i10000").ClearContents
      psynmax = -100;
      psynmaxsh = -100;
      skip = 0; //'this turns off psynthesis

      do //'this loop obtains and stores the entire composite curve
      {
         p = p + 1;
         e = e + einc;
         newtonrhapson(); //'solves for p//'s and e//'s in fingers of chicken foot
                          //'if check >= 400 { //if// GoTo 60: //'NR can//'t find a solution
         if (check > 500) { //if//
            break;//Exit do
         } //End if//
           //'test for total failure
         sum = 0;
         for (k = 1; k <= layers; k++)//k = 1 To layers
         {
            sum = sum + layer[k];
         } //end for k
         if (sum == layers) { //if//
            failspot = "below ground"; //'total failure below ground
            break;//Exit do
         } //End if//
         stem(); //'gets stem and leaf pressures
         leaf();
         if (test == 1)
            break;//Exit do
         compositecurve(); //'stores the entire composite curve
                           //'if skip = 0 { //if//
         leaftemps(); //'gets sun layer leaf temperature from energy balance
         leaftempsshade(); //'gets shade layer leaf temperature
         assimilation(); //'gets sun layer photosynthesis
         assimilationshade(); //'gets shade layer photosynthesis

      } while (!(sum == layers || test == 1 || night == "y" && (dd > 1 && !isNewYear) || check >= 500)); //'loop to complete failure unless it//'s night

      if (chalk > 0) { //if//
         weird = 0; //'done our best
         failspot = "convergence";
         for (z = 1; z <= layers; z++)//z = 1 To layers //'restore layers to functioning if they//'ve been turned off by convergence failure
         {
            if (kminroot[z] != 0)
            {
               layer[z] = 0;
            }
         } //end for z

         goto fortyMarker; //'got as much of the composite curve as is going to happen
      } //End if//
      if (dd == 1 || isNewYear || night == "n") { //if//

         if (check >= 500) { //if// //'try once more
                             //'Stop
            chalk = chalk + 1;

            if (ecritsystem == 0)
            {
               einc = ksatp / 500.0;
               std::cout << "ecritsystem is zero... try resetting to ksatp/500, dd = " << dd << std::endl;
            }

            goto twentyMarker;
         } //End if//

         if (total > 500 || total < 400) { //if//
            einc = ecritsystem / 450.0; //'re-set Einc
            if (ecritsystem == 0)
            {
               einc = ksatp / 500.0;
               std::cout << "ecritsystem is zero... try resetting to ksatp/500, dd = " << dd << std::endl;
            }
            testCount++; // [DEBUG]
            if (testCount > 10)
            {
               testCount = 0;
               goto fortyMarker;
            }
            goto twentyMarker; //'recalculate the composite curve
         } //End if// //'total ok

      } //End if// //'night <>"n" or it//'s not the first round

   fortyMarker:

      bool isNight = true;
      if (night == "n")
         isNight = false;
      if (night == "n" && psynmax > 0 && psynmaxsh > 0 && weird == 0) { //if//
                                                                        //DoEvents //'[HNT] this was required to prevent a hard lock -- this portion of the loop is the most intensive, so let Excel take a "breath" by processing system events to prevent lockup
         canopypressure(); //'returns canopy P and associated output
                           //'if check >= 2000 { //if// GoTo 60:
         if (refilling == "n")
            updatecurves(); //'updates element E(P) curves as required for midday exposure for no refilling
      } //End if// //'night <> "n", psynmax IF

      if (soilred == "y") { //if//
         soilflow(); //'gets vertical soil flow between layers in m3/m2
      }
      else {
         for (z = 0; z <= layers; z++)//z = 0 To layers
         {
            soilredist[z] = 0;
         } //end for z
      } //End if// //'soil red <> y
      if (ground == "y")
         deepflow(); //'gets groundwater flow into bottom layer in m3/m2
                     //'} //End if// //'pet <> y or n
      if (gs_inGrowSeason && sevap == "y") { //if//
         soilevaporation(); //'gets soil evaporation rate
      }
      else {
         soilevap = 0;
      } //End if//

      if (failure == 0 || weird == 1) { //if//
                                        //Debug.Print "DOING A LOOP-8 " & dd

         if (night == "y" || psynmax == 0 || psynmaxsh == 0) { //if// //'set everything to starting point
            k = 0;
            transpiration = eplantl[k]; //'all gas exchange values are for closed stomata
            md = pleaf[k];
            psynact = psyn[k];
            gcmd = gcanw[k]; //'g for water in mmol
            lavpdmd = lavpd[k] * patm;
            cinc = cin[k];
            halt = k;
            transpirationsh = eplantl[k]; //'all gas exchange values are from most recent historical values
            mdsh = pleaf[k];
            psynactsh = psynsh[k];
            gcmdsh = gcanwsh[k]; //'g for water in mmol
            lavpdmdsh = lavpdsh[k] * patm;
            cincsh = cinsh[k];
            haltsh = k; //'halt is index of midday datum
         } //End if// //'night<>y

         dSheet.Cells(rowD + dd, colD + o + dColF_predawn) = pleaf[0]; //'the predawn
                                                                       //'SUN LAYER OUTPUT
         dSheet.Cells(rowD + dd, colD + o + dColF_P) = md; //'the midday
         dSheet.Cells(rowD + dd, colD + o + dColF_E) = transpiration; //'midday transpiration, mmol s-1 m-2 leaf area
         dSheet.Cells(rowD + dd, colD + o + dColF_Gw) = gcmd; //'midday canopy diffusive conductance to water, mmol s-1m-2
         dSheet.Cells(rowD + dd, colD + o + dColF_laVPD) = lavpdmd; //'leaf-to-air vpd
         dSheet.Cells(rowD + dd, colD + o + dColF_leaftemp) = leaftemp[halt]; //'leaf temp
         dSheet.Cells(rowD + dd, colD + o + dColF_ANet) = psynact; //'net A in umol s-1m-2 leaf area
                                                                   //'anet = psynact //'keeping this old setting just in case downstream Anet is needed
                                                                   //'Cells(16 + dd, o + 21) = chalk //'number of NR restarts
         dSheet.Cells(rowD + dd, colD + o + dColF_ci) = cinc * patm * 1000; //'partial pressure of CO2 in Pa
         dSheet.Cells(rowD + dd, colD + o + dColF_PPFD) = qsl; //'umol s-1m-2 photon flux density
                                                               //'SHADE LAYER OUTPUT
         dSheet.Cells(rowD + dd, colD + o + dColF_S_P) = mdsh; //'the midday
         dSheet.Cells(rowD + dd, colD + o + dColF_S_E) = transpirationsh; //'midday transpiration, mmol s-1 m-2 leaf area
         dSheet.Cells(rowD + dd, colD + o + dColF_S_Gw) = gcmdsh; //'midday canopy diffusive conductance to water, mmol s-1m-2
         dSheet.Cells(rowD + dd, colD + o + dColF_S_laVPD) = lavpdmdsh; //'leaf-to-air vpd
         dSheet.Cells(rowD + dd, colD + o + dColF_S_leaftemp) = leaftempsh[haltsh]; //'leaf temp
         dSheet.Cells(rowD + dd, colD + o + dColF_S_Anet) = psynactsh; //'A in umol s-1m-2 leaf area
                                                                       //'anetsh = psynactsh //'see above
                                                                       //'Cells(16 + dd, o + 30) = dpasun //'for debugging purposes
         dSheet.Cells(rowD + dd, colD + o + dColF_S_ci) = cincsh * patm * 1000; //'partial pressure of CO2 in Pa
         dSheet.Cells(rowD + dd, colD + o + dColF_S_PPFD) = qsh; //'umol s-1m-2 photon flux density
                                                                 //'WHOLE TREE OUTPUT
         if (night == "n")
            transpirationtree = laisl / lai * transpiration + laish / lai * transpirationsh; //'weighted mean
         if (night == "y")
            transpirationtree = transpiration;
         if (night == "n")
            atree = laisl / lai * psynact + laish / lai * psynactsh; //'weighted mean
         if (night == "y")
            atree = (psynact + psynactsh) / 2.0; //'simple average at night when there//'s no sun or shade leaves

         dSheet.Cells(rowD + dd, colD + o + dColF_T_E) = transpirationtree; //'weighted mean
         dSheet.Cells(rowD + dd, colD + o + dColF_T_ANet) = atree;
         //'Cells(16 + dd, o + 35) = dpamax //'shade leaf dpa
         //'HYDRAULIC OUTPUT (BASED ON SUN MD)
         dSheet.Cells(rowD + dd, colD + o + dColF_T_pcrit) = pcritsystem;
         dSheet.Cells(rowD + dd, colD + o + dColF_T_Ecrit) = ecritsystem * (1 / laperba) * (1.0 / 3600.0) * 55.4 * 1000; //'ecrit in mmol s-1m-2
         dSheet.Cells(rowD + dd, colD + o + dColF_CP_Pstem) = pstem[halt];
         dSheet.Cells(rowD + dd, colD + o + dColF_CP_Proot) = proot[halt];
         dSheet.Cells(rowD + dd, colD + o + dColF_CP_kstem) = kstem[halt]; //'k stem at midday in kg hr-1m-2MPa-1
         dSheet.Cells(rowD + dd, colD + o + dColF_CP_kleaf) = kleaf[halt]; //'k in leaf at midday

         if (transpiration > 0) { //if//
            kplantold = eplant[halt] / (pleaf[halt] - pleaf[0]);  //'whole plant k at midday in kg hr-1 m-2 basal area...sun value
            dSheet.Cells(rowD + dd, colD + o + dColF_CP_kplant) = kplantold;
         } //End if//
         if (transpiration == 0)
            dSheet.Cells(rowD + dd, colD + o + dColF_CP_kplant) = kplantold; //'use most recent kplant

         if (kplantold < gs_ar_kPlant[gs_yearIndex] || gs_ar_kPlant[gs_yearIndex] == 0)
            gs_ar_kPlant[gs_yearIndex] = kplantold;

         k = o + dColF_CP_kroot1 - 1;//43;
         sum = 0;
         for (z = 1; z <= layers; z++)//z = 1 To layers
         {
            dSheet.Cells(rowD + dd, colD + k + z) = kroot[z][halt]; //'root k at midday, sun fluxes
            sum = sum + kroot[z][halt];
         } //end for z
         dSheet.Cells(rowD + dd, colD + k + 1 + layers) = sum; //'total root k at midday
         if (failure == 0) { //if//
            double tempDouble = 0.0;
            tempDouble = 1 / (1 / kminleaf + 1 / kminstem + 1 / dSheet.Cells(rowD + dd, colD + k + 1 + layers)); //'total xylem k
            dSheet.Cells(rowD + dd, colD + k) = tempDouble; //1 / (1 / kminleaf + 1 / kminstem + 1 / dSheet.Cells(rowD + dd, colD + k + 1 + layers)); //'total xylem k
            if (tempDouble < gs_ar_kXylem[gs_yearIndex] || gs_ar_kXylem[gs_yearIndex] == 0)
               gs_ar_kXylem[gs_yearIndex] = tempDouble;
         } //End if//
         for (z = 1; z <= layers; z++)//z = 1 To layers
         {
            dSheet.Cells(rowD + dd, colD + k + 1 + layers + z) = elayer[z][halt] * (1 / laperba) * (1.0 / 3600.0) * 55.56 * 1000; //'uptake in mmol s-1m-2 leaf area...sun rate
         } //end for z
           //'Cells(16 + dd, k + 2 + 2 * layers) = failspot //'position of failure at critical point

           // TODO since all IO is being handled as double currently, cannot add these failure notes
           // temporarily putting in an obvious number to flag failure
         for (z = 1; z <= 1; z++)//z = 1 To 1
         {
            if (layer[z] == 1)
               dSheet.Cells(rowD + dd, colD + k + 2 + 2 * layers + z) = -1137;// layerfailure[z]; //'layers failed at critical point
         } //end for z

           //Debug.Print "DOING A LOOP-9 " & dd
      } //End if// //'failure IF (basically...failure can//'t happen!)


      if (dd == 1 || isNewYear) { //if// //'NOTE: must be sure that pcritsystem is computed for dd=1!!! (i.e., it//'s not computed at night)
         x = pcritsystem; //'estimate of "permanent wilting point"
         for (z = 1; z <= layers; z++)//z = 1 To layers
         {
            swclimit[z] = swc(x); //'theta/thetasat at critical point
            swclimit[z] = swclimit[z] * thetasat[z]; //'convert to water content
            swclimit[z] = swclimit[z] * depth[z]; //'water content left over in m3/m2 ground area
                                                  //'sumsoil = sumsoil + (fc[z] - swclimit[z]) //'sum is total m3 water per m2 ground withdrawn
         } //end for z

           //Debug.Print "DOING A LOOP-11 " & dd
      } //End if//

        //'now...need to reset layer failure status to midday values (not critical point)
      for (z = 1; z <= layers; z++)//z = 1 To layers
      {
         if (layer[z] == 1) { //if// //'check to see if kminroot[z]=0; otherwise, restore it to function
            if (layerfailure[z] == "root" && kminroot[z] > 0) { //if//
               layer[z] = 0;
               layerfailure[z] = "ok";
            } //End if//
            if (layerfailure[z] == "rhizosphere" && kminroot[z] > 0) { //if//
               layer[z] = 0;
               layerfailure[z] = "ok";
            } //End if//
         } //End if//
      } //end for z

      if (isNewYear)
         isNewYear = false; // always set this

      return -1;

   }

   void readGrowSeasonData()
   {
      long startRow = 2; //we've thrown out the header, but not shifted everything up 1
      long startCol = 1;

      //'assumes table layout:
      //'Year | Start | End, with "Year" label as the named cell above
      long gsC = -1;
      do
      {
         //a[k] = atof(paramCells[rowLR + k][colLR + 3].c_str()); //van genuchten alpha
         gsC = gsC + 1;
         gs_ar_years[gsC] = GSCells[startRow + gsC][startCol]; //gsRange.Offset(gsC + 1, 0).Value2
         gs_ar_starts[gsC] = GSCells[startRow + gsC][startCol + 1]; //gsRange.Offset(gsC + 1, 1).Value2
         gs_ar_ends[gsC] = GSCells[startRow + gsC][startCol + 2]; //gsRange.Offset(gsC + 1, 2).Value2
      } while (!(GSCells[startRow + gsC + 1][startCol] == 0)); // hopefully we get zeros in the blank space...
                                                               //Loop Until IsEmpty(gsRange.Offset(gsC + 1, 0).Value2)
                                                               //growSeasonCount = gsC
   }

   void readIterationSettings()
   {
      if (getValueFromNameStr("i_iter_gwEnable") == "y")
         iter_gwEnable = true;
      else
         iter_gwEnable = false;

      iter_gwInc = getValueFromNameDbl("i_iter_gwInc");
      iter_gwStart = getValueFromNameDbl("i_iter_gwStart");
      iter_gwEnd = getValueFromNameDbl("i_iter_gwEnd");

      if (getValueFromNameStr("i_iter_ffcEnable") == "y")
         iter_ffcEnable = true;
      else
         iter_ffcEnable = false;

      iter_ffcInc = getValueFromNameDbl("i_iter_ffcInc");
      iter_ffcStart = getValueFromNameDbl("i_iter_ffcStart");
      iter_ffcEnd = getValueFromNameDbl("i_iter_ffcEnd");

      if (getValueFromNameStr("i_iter_bagaEnable") == "y")
         iter_bagaEnable = true;
      else
         iter_bagaEnable = false;

      iter_bagaInc = getValueFromNameDbl("i_iter_bagaInc");
      iter_bagaStart = getValueFromNameDbl("i_iter_bagaStart");
      iter_bagaEnd = getValueFromNameDbl("i_iter_bagaEnd");
      iter_bagaRef = getValueFromNameDbl("i_iter_bagaRef");
      iter_bagaCutoff = getValueFromNameDbl("i_iter_bagaCutoff");

      iter_bagaRef = 1.0; //156.269; // 1.0; // TODO TEMP can remove after param sheets updated
      iter_bagaEnd = 500.0; // TODO TEMP " -> for bisection method, allow extreme range

      if (getValueFromNameStr("i_iter_useAreaTable") == "y")
         iter_useAreaTable = true;
      else
         iter_useAreaTable = false;


      if (getValueFromNameStr("i_iter_yearsAsCount") == "y")
         iter_yearsAsCount = true;
      else
         iter_yearsAsCount = false;


      if (getValueFromNameStr("i_iter_runSupplyCurve") == "y")
         iter_runSupplyCurve = true;
      else
         iter_runSupplyCurve = false;

   }

   void modelProgramNewYear()
   {
      // save the iterate-able water system states
      std::string oldGround;
      std::string oldRaining;
      bool oldRainEnabled;
      oldGround = ground;
      oldRaining = raining;
      oldRainEnabled = rainEnabled;

      // clear all the model variables
      cleanModelVars(); //initialize all used variables to avoid any memory effect from previous runs
      readin(); // get all the parameters again

      gs_prevDay = 0;
      gs_inGrowSeason = false;

      // set up the iteration details based on what we saved
      if (iter_runSupplyCurve) {
         ground = "n"; //disable ground water for the first two iterations
         raining = "n"; //disable rain for the first iteration
         rainEnabled = false; //disable rain for the first iteration
         ground = oldGround;
         raining = oldRaining;
         rainEnabled = oldRainEnabled;

         if (iter_bagaEnable) {
            baperga = iter_baga * 0.0001;
            //if we set the BA:GA { also set the LAI
            lai = (laperba * baperga) / treeToPhotoLAI;
            std::cout << "DEBUG: set BA:GA to " << iter_baga << " m2/ha -> " << baperga << "m2/m2. treeLAI/fotoLAI = " << treeToPhotoLAI << " and LAI set to " << lai << std::endl;
         } //End if

      } //End if

      for (k = 0; k <= layers; k++) // To layers //assign source pressures, set layer participation
      {
         layerfailure[k] = "ok";
         layer[k] = 0; //1 if out of function
      } //Next k

        //[HNT] all this is done in C now
      if (true) { //Not useDLL {
         failure = 0; //=1 for system failure at midday...terminates run
         failspot = "no failure";
         componentpcrits(); //gets pcrits for each component
         failspot = "no failure";

         for (k = 1; k <= layers; k++) {//k = 1 To layers //exclude the top layer
            kminroot[k] = ksatr[k];
         } //Next k

         kminstem = ksats;
         kminleaf = ksatl;
         kminplant = ksatp;

         gwflow = 0; //inflow to bottom of root zone
         drainage = 0; //drainage from bottom of root zone
      } //End if


        //if kmaxset = False Or leafpercentset = False { dd = 0 //start with initializing row
        //if kmaxset = True And leafpercentset = True { dd = 1 //skip initializing row
        //dd = 0
        //liveGraphNextUpdate = liveGraphUpdate

        //cleanModelVars(); //initialize all used variables to avoid any memory effect from previous runs
        //readin();
        //Call loadParamsToC
        //Call CPP_setIterationCount(iter_Counter)
   }

   void readSiteAreaValues()
   {
   }

   void saveOutputSheet(std::string filename, std::string sheetType)
   {
      std::cout << "Called saveOutputSheet! Filename " << filename << " type " << sheetType << std::endl;
      // output!
      std::ofstream dataOut;
      dataOut.open(filename + ".csv", std::ios::out);
      if (!dataOut.is_open())
         std::cout << "FAILED to open file " << filename + ".csv" << std::endl;

      dataOut.precision(FIO_PRECISION);

      //int rowCount = -1;
      if (sheetType == "timesteps")
      {
         long maxRow = (gs_yearIndex + 1) * 366 * 24 + 10;
         //if (useGSData)
         //   maxRow = 310000; // 306,600 hours in 35 years (307,440 in 35 leap years) ... leave a buffer

         for (int rowCount = 1; rowCount < maxRow; rowCount++) // excel can only load ~ 1 mil rows
         {
            if (rowCount < rowD)
            {
               for (int hRC = 0; hRC < dataHeaderRow.size(); hRC++)
               {
                  dataOut << dataHeaderRow[hRC] << ",";
               }
            }
            else if (rowCount >= rowD)
            {
               for (int rC = 1; rC < 80; rC++)
               {
                  dataOut << dataCells[rowCount][rC] << ",";
               }
            }
            dataOut << std::endl; // no matter what, write a new line at the end of the row output
         }
         std::cout << "Write complete!" << std::endl;
      }
      else if (sheetType == "summary")
      {
         for (int rowCount = 1; rowCount < MAX_SUMMARY_ROWS; rowCount++) // excel can only load ~ 1 mil rows
         {
            if (rowCount == rowD)
            {
               for (int hRC = 0; hRC < summaryHeaderRow.size(); hRC++)
               {
                  dataOut << summaryHeaderRow[hRC] << ",";
               }
            }
            else if (rowCount > rowD)
            {
               if (finalOutCells[rowCount][dColF_GS_year] > 0.001)
               {
                  for (int rC = 1; rC < MAX_SUMMARY_COLS; rC++)
                  {
                     dataOut << finalOutCells[rowCount][rC] << ",";
                  }
               }
               else
                  break;
            }
            dataOut << std::endl; // no matter what, write a new line at the end of the row output
         }
         std::cout << "Write complete!" << std::endl;
      }
      else if (sheetType == "bagahist")
      {
         // write this one to the root!
         // species/runtype/site/scen/model
         dataOut << stage_OptHistBAGA << "," << std::endl;
         std::cout << "Write complete!" << std::endl;
      }
      else if (sheetType == "bagafut")
      {
         // write this one to the root!
         // species/runtype/site/scen/model
         dataOut << stage_OptFutBAGA << "," << std::endl;
         std::cout << "Write complete!" << std::endl;
      }
   }

   long modelProgramMain(); //program starts here
};

ModelProgram mainProg;
//ModelProgram backup[2];

long ModelProgram::modelProgramMain() //program starts here
{
   //Dim ddOutMod As Long // moved to module global
   memset(finalOutCells, 0, sizeof(finalOutCells)); // clear the final outputs data. Normal data sheets get cleared on each iteration, but this one only per-run

   bool lrSuccess = locateRanges(); //Finds all of the input/output sections across the workbook
                                    // In the C++ version, also loads parameter file and "nametable" for parameters

   if (!lrSuccess)
   {
      std::cout << "Unrecoverable model failure!" << std::endl;
      return 0; // failure, unrecoverable
   }

   readIterationSettings();
   iter_ddOutMod = 0;
   iter_Counter = 0;
   iter_code = 0;

   cleanModelVars();
   initModelVars();
   readin(); //get all global parameters
   if (stage_ID == STAGE_ID_FUT_STRESS_NOACCLIM) // override some if we're doing the odd "no acclimation stress profile"
   {
      std::cout << "Stage " << stage_ID << "; NoAcclim Stress Profile, overriding historical ca " << ca << " -> " << stage_CO2Fut << " and ksatp " << ksatp << " -> " << stage_KmaxFut << std::endl;

   }
   readDataSheet();
   readGSSheet();
   readGrowSeasonData(); //hnt todo cleanup - should just put this in readin?
   if ((iter_useAreaTable)) {
      readSiteAreaValues();
   }

   if (iter_Counter == 0) { //we//re on the first iteration (or we//re not using iterations)

      gs_yearIndex = 0; //multiyear
      gs_prevDay = 0;
      gs_inGrowSeason = false;
   }
   else
   {
      //things to do ONLY if we//re NOT on the first iteration
   } // End If

   memset(gs_ar_input, 0, sizeof(gs_ar_input));
   memset(gs_ar_Anet, 0, sizeof(gs_ar_Anet));
   memset(gs_ar_E, 0, sizeof(gs_ar_E));
   memset(gs_ar_PLCp, 0, sizeof(gs_ar_PLCp));
   memset(gs_ar_PLCx, 0, sizeof(gs_ar_PLCx));
   memset(gs_ar_kPlant, 0, sizeof(gs_ar_kPlant));
   memset(gs_ar_kXylem, 0, sizeof(gs_ar_kXylem));
   memset(gs_ar_ET, 0, sizeof(gs_ar_ET));
   memset(gs_ar_PLC85, 0, sizeof(gs_ar_PLC85));
   memset(gs_ar_PLCSum, 0, sizeof(gs_ar_PLCSum));
   memset(gs_ar_PLCSum_N, 0, sizeof(gs_ar_PLCSum_N));

   memset(gs_ar_kPlantMean, 0, sizeof(gs_ar_kPlantMean));
   memset(gs_ar_kPlantMean_N, 0, sizeof(gs_ar_kPlantMean_N));
   memset(gs_ar_waterInitial, 0, sizeof(gs_ar_waterInitial));
   memset(gs_ar_waterFinal, 0, sizeof(gs_ar_waterFinal));

   memset(gs_ar_waterInitial_GS, 0, sizeof(gs_ar_waterInitial_GS));
   memset(gs_ar_waterFinal_GS, 0, sizeof(gs_ar_waterFinal_GS));
   memset(gs_ar_waterInput_GS, 0, sizeof(gs_ar_waterInput_GS));

   memset(gs_ar_waterInitial_OFF, 0, sizeof(gs_ar_waterInitial_OFF));
   memset(gs_ar_waterFinal_OFF, 0, sizeof(gs_ar_waterFinal_OFF));
   memset(gs_ar_waterInput_OFF, 0, sizeof(gs_ar_waterInput_OFF));

   memset(gs_ar_nrFailConverge, 0, sizeof(gs_ar_nrFailConverge));
   memset(gs_ar_nrFailConverge_Water, 0, sizeof(gs_ar_nrFailConverge_Water));
   memset(gs_ar_nrFailThreshold, 0, sizeof(gs_ar_nrFailThreshold));

   memset(gs_ar_cica, 0, sizeof(gs_ar_cica));
   memset(gs_ar_cica_N, 0, sizeof(gs_ar_cica_N));

   memset(gs_ar_Aci, 0, sizeof(gs_ar_Aci));
   memset(gs_ar_AnetDay, 0, sizeof(gs_ar_AnetDay));

   for (k = 0; k <= layers; k++) // k = 0 To layers //assign source pressures, set layer participation
   {
      layerfailure[k] = "ok";
      layer[k] = 0; //1 if out of function
   } // Next k


   failure = 0; //=1 for system failure at midday...terminates run
   failspot = "no failure";
   componentpcrits(); //gets pcrits for each component
   failspot = "no failure";

   for (k = 1; k <= layers; k++) // k = 1 To layers //exclude the top layer
   {
      kminroot[k] = ksatr[k];
   } // Next k

   kminstem = ksats;
   kminleaf = ksatl;
   kminplant = ksatp;

   gwflow = 0; //inflow to bottom of root zone
   drainage = 0; //drainage from bottom of root zone

   dd = 0;
   long ddMod = 0;
   long successCode = 0;

   do //loop through time steps
   {
      dd = dd + 1;

      successCode = modelTimestepIter(dd);

      if (successCode == 0)
      {
         std::cout << "Unrecoverable model failure!" << std::endl;
         return 0; // failure, unrecoverable
      }
      else if (successCode > 0) // this returns the year if we've incremented it -- not necessary in the full C version (also only supports 1 year right now)
      {
         dd = dd - 1; //we need to repeat this timestep because we bailed early when finding a new year
         gs_yearIndex = successCode;

         // if we're running without growing season limits, we need to record the "end of GS" water content now
         // because we did not complete the previous timestep, back up 1 more to grab a value
         if (!useGSData && gs_yearIndex > 0)
            gs_ar_waterFinal_GS[gs_yearIndex - 1] = dSheet.Cells(rowD + dd - 1, colD + dColF_End_watercontent); // make sure this goes with the previous year

         modelProgramNewYear();
      }
      else // -1 = success, VBA bool convention
      {
         int breakpoint = 1137; // success, in the C version we just continue instead of outputting
                                // do all CSV writing at the end

                                // if we're running without growing season limits, we need to record the "end of GS" water content now
         if (!useGSData && gs_yearIndex > 0)
            gs_ar_waterFinal_GS[gs_yearIndex] = dSheet.Cells(rowD + dd - 1, colD + dColF_End_watercontent); // if this was the end of the set of years, gs_yearIndex will not have been changed so use as-is
      }

      if (dd % 1000 == 0)
         std::cout << "Timestep " << dd << " completed" << std::endl;
   } while (!(dSheet.Cells(rowD + 1 + dd, colD + dColDay) < 0.01)); // loop until the jd value on next row is zero -- it's an integer, but everything is stored in the array as double

                                                                    //Dim gsCount As Long
   long gsCount = 0;
   for (gsCount = 0; gsCount <= gs_yearIndex; gsCount++) // gsCount = 0 To gs_yearIndex
   {
      if (gs_ar_years[gsCount] > 0) { //don//t bother with years that don//t exist
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_year) = gs_ar_years[gsCount]; //gs_ar_Anet(gs_yearIndex) //[HNT] todo improve I don//t like using this hard-coded constant for the size of the C interface array here, maybe check how many data columns there really are in the sheet
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_input) = gs_ar_input[gsCount];
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_Anet) = gs_ar_Anet[gsCount];
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_E) = gs_ar_E[gsCount];
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_PLCp) = gs_ar_PLCp[gsCount];
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_PLCx) = gs_ar_PLCx[gsCount];
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_kPlant) = gs_ar_kPlant[gsCount];
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_kXylem) = gs_ar_kXylem[gsCount];
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET) = gs_ar_ET[gsCount];
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 1) = rainEnabled;
         if (ground == "y")
            dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 2) = 1.0;
         else
            dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 2) = 0.0;
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 3) = ffc;
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 4) = grounddistance;
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 5) = gs_ar_PLC85[gsCount];
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 6) = baperga / 0.0001;
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 7) = laperba; // need this for the final calcs
                                                                                                                          //dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 7) = gs_ar_kPlantMean[gsCount] / gs_ar_kPlantMean_N[gsCount];
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 8) = gs_ar_PLCSum[gsCount] / gs_ar_PLCSum_N[gsCount];
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 9) = gs_ar_waterFinal[gsCount] - gs_ar_waterInitial[gsCount];
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 10) = gs_ar_waterInitial[gsCount]; // we want to know what the initial was too, in case it's not FC

         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 11) = gs_ar_waterInitial_OFF[gsCount]; // initial content for preceding off-season
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 12) = gs_ar_waterInput_OFF[gsCount]; // input for preceding off-season
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 13) = gs_ar_waterFinal_OFF[gsCount]; // final content for preceding off-season

         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 14) = gs_ar_waterInitial_GS[gsCount]; // initial content for growing season
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 15) = gs_ar_waterInput_GS[gsCount]; // input for growing season
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 16) = gs_ar_waterFinal_GS[gsCount]; // final content for growing season

         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 17) = gs_ar_nrFailConverge[gsCount]; // number of convergence failures
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 18) = gs_ar_nrFailConverge_Water[gsCount] / gs_ar_nrFailConverge[gsCount]; // avg water content during convergence failure
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 19) = gs_ar_nrFailConverge_WaterMax[gsCount]; // MAX water content during convergence failure
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 20) = gs_ar_nrFailThreshold[gsCount]; // MAX water content during convergence failure

         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 21) = gs_ar_cica[gsCount] / gs_ar_cica_N[gsCount];

         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 23) = (gs_ar_Aci[gsCount] / gs_ar_AnetDay[gsCount]) * patm * 1000.0;
         dSheet.fCells(rowD + iter_Counter + 1 + gsCount /* *40 */, colD /*+ gsCount * 16*/ + dColF_GS_ET + 24) = (gs_ar_Aci[gsCount] / gs_ar_AnetDay[gsCount]) / ca;

         /*memcpy(&testProg, this, sizeof(ModelProgram));
         std::cout << "TEST! My REAL ci/ca " << gs_ar_cica[gsCount] / gs_ar_cica_N[gsCount] << std::endl;
         std::cout << "TEST! My COPIED ci/ca " << testProg.gs_ar_cica[gsCount] / testProg.gs_ar_cica_N[gsCount] << std::endl;*/
      }
      else
      {
         break; //quit the loop if we reach the end of the growing seasons list early somehow
      } // End If
   } // Next gsCount
   saveOutputSheet("./" + stageNames[stage_ID] + "_OUTPUT_timesteps", "timesteps");
   saveOutputSheet("./" + stageNames[stage_ID] + "_OUTPUT_summary", "summary");

   return 1;
}

int main()
{
   // seed the random number generator with something crazy
   srand((unsigned)(time(0) * time(0)));
   // set the cout decimal precision
   std::cout.precision(12);

   stageNames[STAGE_ID_NONE] = "standard"; //stage IDs and names are related to various modes the model was designed to run in for specific projects -- that code has been removed, but a few vestiges remain

   long result = 0;
   std::cout << "Model Startup." << std::endl;

   // to do a normal run that's only based on local folder parameter sheet settings, set the stage_ID to zero
   mainProg.stage_ID = STAGE_ID_NONE;
   result = mainProg.modelProgramMain(); // for returning failure error codes... not really used in this version

   if (!result)
   {
      std::cout << "Model Failure! Stage " << mainProg.stage_ID << std::endl;
      return 0;
   }
   else
   {
      std::cout << "Model Success! Stage " << mainProg.stage_ID << std::endl;
   }
   return 1;
}