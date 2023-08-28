#include <cstdio>
#include <string>
#include <cstring>
#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <climits>

#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Io/FileTools.h>


using namespace bpp;
using namespace std;


struct dfem_data
{
  int nb_gene;
  bool folded;
  int nb_class;
  double* SFS_S;
  double* SFS_N;
  double dN, dNa, dNna, dS;
  double LN, LS;

  double true_dnds;
  double true_omegaA;
  double true_omegaNA;
  double true_alpha;
};


/* usage */

void usage()
{
  cerr << "usage: slim_to_dfem SLIM_outfput_file DFEM_input_file folded|unfolded" << endl << endl;
  exit(0);
}


/* read_slim */

struct dfem_data read_slim(vector<string> slimv, bool folded)
{
  struct dfem_data dd;
  int burnin;

  // check mutation types
  int unsigned i = 0;
  while (slimv[i] != "#MUTATION TYPES")
  {
    i++;
    if (i == slimv.size())
    {
      cerr << "missing info on #MUTATION TYPES in slim file" << endl; exit(0);
    }
  }
  if (i + 3 >= slimv.size())
  {
    cerr << "3 mutation types expected in slim file" << endl; exit(0);
  }
  if (slimv[i + 1].find("initializeMutationType(1") == string::npos)
  {
    cerr << "missing mutation type m1 in slim file" << endl; exit(0);
  }
  if (slimv[i + 2].find("initializeMutationType(2") == string::npos)
  {
    cerr << "missing mutation type m2 in slim file" << endl; exit(0);
  }
  if (slimv[i + 3].find("initializeMutationType(3") == string::npos)
  {
    cerr << "missing mutation type m3 in slim file" << endl; exit(0);
  }
  cerr << "# m1 mutations are considered neutral" << endl;
  cerr << "# m2 mutations are considered selected (deleterious)" << endl;
  cerr << "# m3 mutations are considered selected (advantageous)" << endl;

  // nb of sites
  while (slimv[i] != "#GENOMIC ELEMENT TYPES")
  {
    i++;
    if (i == slimv.size())
    {
      cerr << "missing info on #GENOMIC ELEMENT TYPES in slim file" << endl; exit(0);
    }
  }
  i++;
  StringTokenizer* tok = new StringTokenizer(slimv[i], " \n");
  std::deque< std::string > gen = tok->getTokens();
  delete(tok);
  if (gen.size() != 3)
  {
    cerr << "3 tokens expected in GENOMIC ELEMENT TYPES line" << endl; exit(0);
  }
  if (gen[1] != "m1")
  {
    cerr << "2nd token in GENOMIC ELEMENT TYPES line should be m1" << endl; exit(0);
  }
  double S_site_prop = TextTools::toDouble(gen[2]);
  i++;

  std::deque< std::string > chr;
  int total_lg = 0;
  while (slimv[i] != "#RECOMBINATION RATE")
  {
    i++;
    if (i == slimv.size())
    {
      cerr << "missing info on #GENOMIC ELEMENT TYPES in slim file" << endl; exit(0);
    }
    if (slimv[i][0] == '#')
      continue;
    tok = new StringTokenizer(slimv[i], ", );\n");
    chr = tok->getTokens();
    delete(tok);
    if (chr[0] != "initializeGenomicElement(g1")
    {
      cerr << "a single GENOMIC ELEMENT TYPE (g1) expected" << endl; exit(0);
    }
    total_lg += TextTools::toInt(chr[2]) - TextTools::toInt(chr[1]);
  }

  dd.LS = S_site_prop * total_lg;
  dd.LN = total_lg - dd.LS;

  // nb of genes & classes
  dd.folded = folded;
  i = 0;
  while (slimv[i].find("#OUT:") == string::npos)
  {
    i++;
    if (i == slimv.size())
    {
      cerr << "missing output in slim file" << endl; exit(0);
    }
  }
  tok = new StringTokenizer(slimv[i], " \n");
  std::deque< std::string > out_l1 = tok->getTokens();
  delete(tok);
  if (out_l1.size() != 7)
  {
    cerr << "6 tokens expected in first output #OUT: line" << endl; exit(0);
  }
  burnin = TextTools::toInt(out_l1[1]);
  dd.nb_gene = TextTools::toInt(out_l1[5]);
  if (!folded)
    dd.nb_class = dd.nb_gene - 1;
  else
    dd.nb_class = dd.nb_gene / 2;

  // allocate
  dd.SFS_S = (double*)calloc(dd.nb_class, sizeof(double));
  dd.SFS_N = (double*)calloc(dd.nb_class, sizeof(double));

  // count polymorphisms
  ++i;
  while (slimv[i].find("#OUT:") == string::npos)
  {
    i++;
    if (i == slimv.size())
    {
      cerr << "missing output in slim file" << endl; exit(0);
    }
  } 
  bool count = true;
  std::deque< std::string > mut;
  while (slimv[i][slimv[i].size() - 1] != 'F')
  {
    i++;
    if (slimv[i] == "Mutations:")
    {
      count = true; continue;
    }
    if (slimv[i] == "Genomes:")
    {
      count = false; continue;
    }
    if (slimv[i].find("#OUT:") != string::npos)
      continue;
    if (count)
    {
      tok = new StringTokenizer(slimv[i], " \n");
      mut = tok->getTokens();
      delete(tok);
      if (mut.size() != 9)
      {
        cerr << "9 tokens expected in mutation output lines" << endl << slimv[i] << endl; exit(0);
      }
      int freq = TextTools::toInt(mut[8]);
      if (freq > dd.nb_gene)
      {
        cerr << "unexpected mutation frequency: " << freq << endl << slimv[i] << endl; exit(0);
      }
      if (freq == dd.nb_gene)
        continue;
      int cla;
      if (!folded)
        cla = freq - 1;
      if (folded)
      {
        if (freq <= dd.nb_class)
          cla = freq - 1;
        else
          cla = dd.nb_gene - freq - 1;
      }
      if (mut[2] == "m1")
        dd.SFS_S[cla]++;
      else if (mut[2] == "m2" || mut[2] == "m3")
        dd.SFS_N[cla]++;
    }
  }

  // count substitutions
// burnin: échantillonner à 10*N pour que ce soit le premier #OUT et corresponde à la période de burnin au cours delaquelle on ne compte pas les substitutions.
  double d1, d2, d3;
  d1 = d2 = d3 = 0;
  std::deque< std::string > subst;
  while (i < slimv.size() - 1)
  {
    i++;
    if (slimv[i] == "Mutations:")
      continue;
    if (slimv[i] == "")
      continue;
    tok = new StringTokenizer(slimv[i], " \n");
    subst = tok->getTokens();
    delete(tok);
    if (subst.size() != 9)
    {
      cerr << "9 tokens expected in substitution output lines" << endl << slimv[i] << endl; exit(0);
    }
    int gen = TextTools::toInt(subst[8]);
    if (gen < burnin)
      continue;
    if (subst[2] == "m1")
      d1++;
    else if (subst[2] == "m2")
      d2++;
    else if (subst[2] == "m3")
      d3++;
  }
  dd.dNa = d3;
  dd.dNna = d2;
  dd.dN = d2 + d3;
  dd.dS = d1;

  dd.true_dnds = (dd.dN / dd.LN) / (dd.dS / dd.LS);
  dd.true_omegaA = (dd.dNa / dd.LN) / (dd.dS / dd.LS);
  dd.true_omegaNA = (dd.dNna / dd.LN) / (dd.dS / dd.LS);
  dd.true_alpha = d3 / dd.dN;

  cout << "# d1:" << d1 << " d2:" << d2 << " d3:" << d3 << endl;
  cout << "# LN:" << dd.LN << " LS:" << dd.LS << endl;

  return dd;
}


/*************************************/
/************    MAIN    *************/
/*************************************/

int main(int argc, char** argv)
{
  string slimfile_name, dfefile_name;
  std::ifstream* slimfile;
  std::ofstream* fout;
  vector<string> slimv;
  bool folded;
  struct dfem_data dd;

  // read arguments
  if (argc != 4)
    usage();
  slimfile_name = argv[1];
  dfefile_name = argv[2];
  if (strcmp(argv[3], "folded") == 0)
    folded = true;
  else
    folded = false;

  if (fopen(slimfile_name.c_str(), "r") == NULL)
  {
    cerr << "can't open slim file" << endl; exit(0);
  }

  // open files
  slimfile = new std::ifstream(slimfile_name.c_str(), ios::in);
  fout = new std::ofstream(dfefile_name.c_str(), ios::out);

  // put content into vectors of string
  slimv = FileTools::putStreamIntoVectorOfStrings(*slimfile);

  cerr << "# slim file has " << slimv.size() << " lines" << endl;

  // read slim file
  dd = read_slim(slimv, folded);

  // write dfem file
  *fout << "from SLIM output file " << slimfile_name << endl;
  if (!folded)
    *fout << "#unfolded" << endl;
  *fout << "all" << "\t" << dd.nb_gene << "\t";
  *fout << dd.LN << "\t";
  for (int ii = 0; ii < dd.nb_class; ii++)
  {
    *fout << dd.SFS_N[ii] << "\t";
  }
  *fout << dd.LS << "\t";
  for (int ii = 0; ii < dd.nb_class; ii++)
  {
    *fout << dd.SFS_S[ii] << "\t";
  }
  *fout << dd.LN << "\t" << dd.dN << "\t" << dd.LS << "\t" << dd.dS << endl;

  // summary
  double totN, totS;
  totN = totS = 0;
  for (size_t ii = 0; ii < dd.nb_class; ++ii)
  {
    totN += dd.SFS_N[ii];
    totS += dd.SFS_S[ii];
  }
  cout << endl;
  cout << "# found:" << endl;
  cout << "#  " << totN << " selected SNPs" << endl;
  cout << "#  " << totS << " neutral SNPs" << endl;
  cout << "#  " << dd.dN << " selected substitutions" << endl;
  cout << "#  " << dd.dS << " neutral substitutions" << endl;
  cout << "#  " << dd.LN << " selected sites" << endl;
  cout << "#  " << dd.LS << " neutral sites" << endl;
  cout << endl;
  cout << "# true (=simulated) dN/dS: " << dd.true_dnds << endl;
  cout << "# true (=simulated) alpha: " << dd.true_alpha << endl;
  cout << "# true (=simulated) omegaA: " << dd.true_omegaA << endl;
  cout << "# true (=simulated) omegaNA: " << dd.true_omegaNA << endl;
  cout << endl;
  cout << "# Polymorphism and divergence data written to file: " << dfefile_name << endl << endl;
  cout << "SelectedSNPs,NeutralSNPs,SelectedSubstitutions,NeutralSubstitutions,SelectedSites,NeutralSites,dNdS,omegaA,omegaNA,alpha" << endl;
  cout << totN << "," << totS << "," << dd.dN << "," << dd.dS << "," << dd.LN << "," << dd.LS << "," << dd.true_dnds << "," << dd.true_omegaA << "," << dd.true_omegaNA << "," << dd.true_alpha << endl;
}
