#ifndef CGVSIMULATION_H_
#define CGVSIMULATION_H_

#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <time.h>
#include <algorithm>
#include <string>
#include <math.h>
#include <set>
#include <map>
#include <vector>
using namespace std;

//test code
extern int total_len1;
extern int total_len2;
//test code
extern double ins_len_ratio[20];
extern double del_len_ratio[20];
const int flank=10;
const int kLargeflank = 1000;
extern bool VAR_type[5];
//const double indel_len_ratio[20]={0.43,0.52,0.87,0.89,0.90,0.94,0.947,0.952,0.972,0.975,0.977,0.987,0.988,0.9888,0.9948,0.9957,0.9997,0.9999,0.99996};
const string snp_set[16]={"AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"};

const double A_ratio[16] = {0/90000.0, 6428/90000.0, 19284/90000.0, 25712/90000.0, 32140/90000.0, 38568/90000.0, 38576/90000.0, 38580/90000.0, 51427/90000.0, 51435/90000.0, 77139/90000.0, 77147/90000.0, 83569/90000.0, 83573/90000.0, 83581/90000.0, 1.0};
const double T_ratio[16] = {6428/90000.0, 6436/90000.0, 6440/90000.0, 12859/90000.0, 12867/90000.0, 38571/90000.0, 38579/90000.0, 51429/90000.0, 51433/90000.0, 51441/90000.0, 57860/90000.0, 64288/90000.0, 70716/90000.0, 83572/90000.0, 90000/90000.0, 1.0};
const double C_ratio[16] = {6428/90000.0, 12856/90000.0, 12860/90000.0, 12868/90000.0, 19287/90000.0, 19287/90000.0, 25715/90000.0, 38571/90000.0, 38575/90000.0, 45000/90000.0, 51428/90000.0, 51436/90000.0, 51444/90000.0, 64288/90000.0, 64296/90000.0, 1.0};
const double G_ratio[16] = {25710/90000.0, 25718/90000.0, 38568/90000.0, 38576/90000.0, 38584/90000.0, 45000/90000.0, 51428/90000.0, 51432/90000.0, 64285/90000.0, 70713/90000.0, 70713/90000.0, 77141/90000.0, 77149/90000.0, 77153/90000.0, 83572/90000.0, 1.0};

//const double blood_snp_variation_ratio=static_cast <double >(50)/57;
//const double blood_indel_variation_ratio=static_cast <double >(6)/57;
//const double blood_sv_variation_ratio=static_cast <double >(2)/1101;

//const double para_snp_variation_ratio=static_cast <double >(49)/61;//60
//const double para_indel_variation_ratio=static_cast <double >(10)/61;//20
//const double para_sv_variation_ratio=static_cast <double >(1)/61;//1

//const double cancer_snp_variation_ratio=static_cast <double >(49)/61;//600
//const double cancer_indel_variation_ratio=static_cast <double >(10)/61;//354
//const double cancer_sv_variation_ratio=static_cast <double >(1)/61;//6

//long
const int kInversionLenL = 1000;
const int kInversionLenU = 200000;
const int kIndelLenL = 1000;
const int kIndelLenU = 4000000;
const int kCopyTimes = 2;

const int kTranslocationDenominator = 4;
const int kTranslocationNumerator = 1;

const double kVirusIntegrationRationInInsertion = 0.1;

//for deletion
//const int kDelDenominator = 4;
//const int kDel

struct var_info
{
	//string variation_type; "ins","del","snp"
	string var_type;
	char ref;
	string seq1;
	string seq2;
	int length1;
	int length2;
	char period;
};

struct sv_var_info
{
    //variation type
    string var_type;
    int allele;
    int length;
	char period;
};

struct insertion_info
{
	string chr;
	map<int, sv_var_info>::iterator pos2sv_var_info_it;
};

struct MVSCOption
{
	string hg_ref;
	string virus_ref;
	string output_dir;
	string simu_str;//five number to Initialization 
	string chr;
	string narea_file;
	string blood_ref;
	string para_ref;
	string cancer_ref;
	/*
	string blood_info;
	string para_info;
	string cancer_info;
	*/
	string blood_snp_info;
	string blood_indel_info;
	string blood_sv_info;
	string blood_inseq_info;
	string blood_cnv_info;
	string para_snp_info;
	string para_indel_info;
	string para_sv_info;
	string para_inseq_info;
	string para_cnv_info;
	string cancer_snp_info;
	string cancer_indel_info;
	string cancer_sv_info;
	string cancer_inseq_info;
	string cancer_cnv_info;

	double blood_variation_ratio;
	double insertion_rate;
	double inversion_rate;
	double deletion_rate;
	double snp_variation_ratio;
	double indel_variation_ratio;
	double sv_variation_ratio;
	double cancer_sv_variation_ratio;
	double para_variation_ratio;
	double cancer_variation_ratio;
	double retain_paraTocancer_ratio;
	//int times;
	MVSCOption()
	{//Initialization
		chr.assign("chr17");


		blood_variation_ratio=0.0006;
		snp_variation_ratio=50/57.0;
		indel_variation_ratio=6/57.0;
		sv_variation_ratio=2/1101.0;
		cancer_sv_variation_ratio=10*sv_variation_ratio;
		para_variation_ratio=0.00001;
		cancer_variation_ratio=0.00004;
		insertion_rate=0.3;
		deletion_rate=0.45;
		inversion_rate=0.25;
		retain_paraTocancer_ratio=0.8;
	}
	void Compute()
	{
		blood_ref.assign(output_dir).append("/germline.fa");
		blood_snp_info.assign(output_dir).append("/germline.snp.info");
		blood_indel_info.assign(output_dir).append("/germline.indel.info");
		blood_sv_info.assign(output_dir).append("/germline.sv.info");
		blood_inseq_info.assign(output_dir).append("/germline.insertseq.fa");
		blood_cnv_info.assign(output_dir).append("/germline.cnv.info");
		//blood_info.assign(output_dir).append("/blood.info");

		para_ref.assign(output_dir).append("/para.fa");
		para_snp_info.assign(output_dir).append("/para.snp.info");
		para_indel_info.assign(output_dir).append("/para.indel.info");
		para_sv_info.assign(output_dir).append("/para.sv.info");
		para_inseq_info.assign(output_dir).append("/para.insertseq.fa");
		para_cnv_info.assign(output_dir).append("/para.cnv.info");
	//	para_info.assign(output_dir).append("/para.info");
		cancer_ref.assign(output_dir).append("/cancer.fa");
		cancer_snp_info.assign(output_dir).append("/cancer.snp.info");
		cancer_indel_info.assign(output_dir).append("/cancer.indel.info");
		cancer_sv_info.assign(output_dir).append("/cancer.sv.info");
		cancer_inseq_info.assign(output_dir).append("/cancer.insertseq.fa");
		cancer_cnv_info.assign(output_dir).append("/cancer.cnv.info");
	//	cancer_info.assign(output_dir).append("/cancer.info");
	}
};
struct variation_ratio
{
	double b_snp_variation_ratio;
	double b_indel_variation_ratio;
	double b_ins_variation_ratio;
	double b_del_variation_ratio;
	double b_inv_variation_ratio;
	double p_snp_variation_ratio;
	double p_indel_variation_ratio;
	double p_ins_variation_ratio;
	double p_del_variation_ratio;
	double p_inv_variation_ratio;
	double c_snp_variation_ratio;
	double c_indel_variation_ratio;
	double c_ins_variation_ratio;
	double c_del_variation_ratio;
	double c_inv_variation_ratio;
	void compute(const MVSCOption& mvsc)
	{
		b_indel_variation_ratio=mvsc.blood_variation_ratio*mvsc.indel_variation_ratio;
		b_snp_variation_ratio=mvsc.blood_variation_ratio*mvsc.snp_variation_ratio;
		b_ins_variation_ratio=mvsc.blood_variation_ratio*mvsc.sv_variation_ratio*mvsc.insertion_rate;
		b_del_variation_ratio=mvsc.blood_variation_ratio*mvsc.sv_variation_ratio*mvsc.deletion_rate;
		b_inv_variation_ratio=mvsc.blood_variation_ratio*mvsc.sv_variation_ratio*mvsc.inversion_rate;
		p_indel_variation_ratio=mvsc.para_variation_ratio*mvsc.indel_variation_ratio;
		p_snp_variation_ratio=mvsc.para_variation_ratio*mvsc.snp_variation_ratio;
		p_ins_variation_ratio=mvsc.para_variation_ratio*mvsc.sv_variation_ratio*mvsc.insertion_rate;
		p_del_variation_ratio=mvsc.para_variation_ratio*mvsc.sv_variation_ratio*mvsc.deletion_rate;
		p_inv_variation_ratio=mvsc.para_variation_ratio*mvsc.sv_variation_ratio*mvsc.inversion_rate;
		c_indel_variation_ratio=mvsc.cancer_variation_ratio*mvsc.indel_variation_ratio;
		c_snp_variation_ratio=mvsc.cancer_variation_ratio*mvsc.snp_variation_ratio;
		c_ins_variation_ratio=mvsc.cancer_variation_ratio*mvsc.cancer_sv_variation_ratio*mvsc.insertion_rate;
		c_del_variation_ratio=mvsc.cancer_variation_ratio*mvsc.cancer_sv_variation_ratio*mvsc.deletion_rate;
		c_inv_variation_ratio=mvsc.cancer_variation_ratio*mvsc.cancer_sv_variation_ratio*mvsc.inversion_rate;

	}
};
//record the variation infomation
//static set<int> point_var;
extern map<int,  var_info> pos2var_info;
extern map<int, sv_var_info> pos2sv_info;
extern map<string, map<int, int> > chr2nbg2nend;
extern map<pair<int, int>, pair<int,char> > pos2allele;
extern map<pair<int, int>, map<pair<int, int>, pair<int,char> > > pos2cppos_allele;

//function
void OpenFiles(MVSCOption mvscoption, ofstream &FOUT_B, ofstream &FOUTSNP_B, ofstream &FOUTINDEL_B, ofstream &FOUTSV_B, ofstream &FOUTSEQ_B, ofstream &FOUTCNV_B, ofstream &FOUT_P, ofstream &FOUTSNP_P, ofstream &FOUTINDEL_P, ofstream &FOUTSV_P, ofstream &FOUTSEQ_P, ofstream &FOUTCNV_P, ofstream &FOUT_C, ofstream &FOUTSNP_C, ofstream &FOUTINDEL_C, ofstream &FOUTSV_C, ofstream &FOUTSEQ_C, ofstream &FOUTCNV_C);
void Simulation_start(map<string, string> &chr2seq, vector<string> &vec_chr, map<string, string> &virus2seq, vector<string> &virus_vec, map<string, map<int, var_info> >& chr2pos2var_info, map<string, map<int, sv_var_info> >& chr2pos2sv_info, map<string, map<int, sv_var_info> >& virus2pos2sv_info, map<string, map<int, insertion_info> > &chr2pos2insertion_info, map<string, map<int, int> >& chr2nbg2nend, const variation_ratio& var, char level);
void Small_variation(const string &, const string& ,const variation_ratio& , map<string, map<int, var_info> >& ,char );
void Long_variation(const string &chr, const string& hg_seq, const variation_ratio& var, map<string, string> &virus2seq, map<string, map<int, sv_var_info> >& chr2pos2sv_info, map<string, map<int, insertion_info> > &chr2pos2insertion_info, map<string, map<int, int> >& chr2nbg2nend, char level);
void Record_deletion(double ,map<int,var_info>& ,map<int, sv_var_info>& ,map<pair<int, int>,  pair<int, char> >& ,map<pair<int, int>, map<pair<int, int>,  pair<int, char> > >& );
void small_var_deletion(double ,map<int,var_info>& );
void long_var_deletion(double ,map<int, sv_var_info>& ,map<pair<int, int>,  pair<int, char> >& ,map<pair<int, int>, map<pair<int, int>,  pair<int, char> > >& );
void Output(ofstream &FOUT, ofstream &FOUTSNP, ofstream &FOUTINDEL, ofstream &FOUTSV, ofstream &FOUTSEQ, ofstream &FOUTCNV, string chr, string outseq1,string outseq2,map<int, var_info> pos2var_info,map<int, sv_var_info> pos2sv_info,map<pair<int, int>, pair<int, char> > pos2allele,map<pair<int, int>, map<pair<int, int>, pair<int, char> > > pos2cppos_allele);
void Output_info(map<int, var_info>& ,map<int, sv_var_info>& ,map<pair<int, int>, pair<int, char> >& ,map<pair<int, int>, map<pair<int, int>, pair<int, char> > >& ,string );
void Output_small_info(map<int, var_info>& pos2var_info, ofstream &FOUTSNP, ofstream &FOUTINDEL, string chr);
void Output_large_info(map<int, sv_var_info>& pos2sv_info, map<pair<int, int>, pair<int, char> >& pos2allele, map<pair<int, int>, map<pair<int, int>, pair<int, char> > >& pos2cppos_allele, ofstream &FOUTSV, ofstream &FOUTSEQ, ofstream &FOUTCNV, string chr);
void OutputDelInvInfo(ofstream &FOUTSV, map<int, sv_var_info>& pos2sv_info, string chr);
void OutputSeq(ofstream &FOUT, string seq);
void OutputFasta(ofstream &FOUT, string chr, string seq);
void OutputInsTranInfo(ofstream &FOUTSV, map<int, insertion_info>& pos2insertion_info, string chr);
//small varation
void snp_simulation(const string &, const string& , double , map<string, map<int, var_info> >& ,char );
void indel_simulation(const string &, const string& , double , map<string, map<int,var_info> >& ,char );
//int get_insert_length();
int get_ins_length();
int get_del_length();
int get_ref_position(int ,int);
string SNP_var(const char &);
//long varation
//void inversion_simulation(const string &, const string& , double , map<string, map<int, sv_var_info> >& , map<int, int>& , char );
//void deletion_simulation(const string &, const string& , double , map<string, map<int, sv_var_info> >& , map<int, int>& , map<pair<int, int>, pair<int, char> >& , char );
//void deletion_simulation(const string &chr, const string& seq, double del_variation_ratio, map<string, string> &id2seq, vector<string> &id_vec, map<string, map<int, insertion_info> > &chr2pos2insertion_info,  map<string, map<int, sv_var_info> >& chr2pos2sv_info, map<string, map<int, int> >& chr2nbg2nend, map<pair<int, int>, pair<int, char> >& pos2allele, char level);
//void insertion_simulation(const string &, const string& , double , map<string, string> & , map<string, map<int, sv_var_info> >& , map<string, map<int, int> >& , map<pair<int, int>, map<pair<int, int>, pair<int, char> > >& , char );
void insertion_simulation(map<string, string> &id2seq, vector<string> &id_vec, map<string, map<int, insertion_info> > &chr2pos2insertion_info, map<string, map<int, sv_var_info> >& chr2pos2sv_info, map<string, map<int, int> >& chr2nbg2nend, string chr_of_insert_seq, map<int, sv_var_info>::iterator &sv_var_info_it);

void VirusIntegrationSimulation(map<string, string> &chr2seq, map<string, string> &virus2seq, vector<string> &virus_vec, map<string, map<int, insertion_info> > &chr2pos2insertion_info, map<string, map<int, sv_var_info> >& chr2pos2sv_info, map<string, map<int, sv_var_info> >& virus2pos2sv_info, map<string, map<int, int> >& chr2nbg2nend, double variation_ratio, char level);
//for large deletion, large copy , large inversion
void del_copy_inv_simulation(const string &chr, const string& seq, double variation_ratio, map<string, map<int, sv_var_info> >& chr2pos2sv_info, map<string, map<int, insertion_info> > &chr2pos2insertion_info, map<string, map<int, int> >& chr2nbg2nend, string type, char level);
void variation_integration(const string &, map<int, var_info> &, map<int, sv_var_info> &, string &, string &);
void VariationIntegration(ofstream &FASTA, ofstream &FOUTSNP, ofstream &FOUTINDEL, ofstream &FOUTSV, map<string, string> &chr2seq, map<string, string> &virus2seq, vector<string> &vec_chr, vector<string> &vec_virus, map<string, map<int, var_info> > &chr2pos2var_info, map<string, map<int, insertion_info> > &chr2pos2insertion_info, map<string, map<int, sv_var_info> >& chr2pos2sv_info, map<string, map<int, sv_var_info> >& virus2pos2sv_info, map<string, map<int, int> >& chr2nbg2nend, char level);

void sv_external_integration(const string& , map<int, var_info>& , map<int, var_info>::iterator& , int& , int , string& , string& );
//void sv_internal_integration(const string & , map<int, var_info>& , map<int, var_info>::iterator& , map<int, sv_var_info>::iterator& , int& , int , string& , string& );
void sv_internal_integration(const string& seq, map<int, var_info>& pos2var_info, map<int, var_info>::iterator& snp_indel_map_it, map<int, sv_var_info>::iterator& sv_map_it, int& last_pos, string& outseq1, string& outseq2);
void ldel_internal_integration(const string& , map<int, var_info>& , map<int, var_info>::iterator& , int& , int , string& , int );

void DeleteSmallDelSpanSpecificPos(map<int, var_info> &pos2var_info, int pos);

//void InsertionIntegration(const string& seq, map<int, var_info>& pos2var_info, map<int, insertion_info>::iterator &pos2insertion_info_it, string& outseq1, string& outseq2);
//void InsertionIntegration(map<string, string>& chr2seq, map<string, map<int, var_info> >& chr2pos2var_info, map<int, insertion_info>::iterator &pos2insertion_info_it, string& outseq1, string& outseq2);
//void InsertionIntegration(map<string, string>& chr2seq, map<string, map<int, var_info> >& chr2pos2var_info, map<string, map<int, var_info> >& virus2pos2var_info, map<int, insertion_info>::iterator &pos2insertion_info_it, string& outseq1, string& outseq2);
void InsertionIntegration(map<string, string>& chr2seq, map<string, string> &virus2seq, map<string, map<int, var_info> >& chr2pos2var_info, map<string, map<int, sv_var_info> >& virus2pos2sv_var_info, map<int, insertion_info>::iterator &pos2insertion_info_it, string& outseq1, string& outseq2);

void ins_cp_simulation(const string &, const string& , map<int, sv_var_info>& , map<string, map<int, int> >& , map<pair<int, int>,map<pair<int, int>, pair<int, char> > >& , char );
void inverse_seq(string& );
bool isInNArea(int , int , map<int, int>& );
bool isInNArea(string, int , int , map<string, map<int, int> >& );
bool isInNAreadForIns(int , map<int, int>& );
bool isInNAreadForIns(string, int , map<string, map<int, int> >& );
//bool isOverlap(int pos, int , map<int, sv_var_info>& );
//bool isInVariationPos(int , map<int, sv_var_info>& );
void display_usage();
double uni_random(double minValue,double maxValue);
double beta_random();
bool isCanInsert(int pos, int len, map<int, var_info> &pos2var_info);
bool isCanInsert(int pos, int len, map<int, sv_var_info> &pos2sv_info);
bool isCanInsert(string chr, int pos, int len, map<string, map<int, sv_var_info> > &chr2pos2sv_info);
bool isLocked(int pos, int len, map<int, insertion_info> &pos2insertion_info);
bool isSmallAndSmallOverlap(int pos, int len, map<int, var_info>::iterator bound_map_it);
bool isLargeAndLargeOverlap(int pos, int len, map<int, sv_var_info>::iterator bound_map_it);
bool is2AreaOverlap(int alow, int aup, int blow, int bup);

int ReadVirusSeq(string file, map<string, string> &virus2seq, vector<string> &vec_virus);
string GetStringBeforeSpace(string str);
				

#endif
