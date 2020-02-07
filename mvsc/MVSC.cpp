#include "cgvsimulation.h"
#include "cgvsimulation.cpp"
#include "SNP_INDELsimulation.cpp"
#include<getopt.h>
#include <ctype.h>
	
bool VAR_type[5] = {1, 1, 1, 1, 1};
map<string, map<int,  var_info> > chr2pos2var_info;
map<string, map<int,  sv_var_info> > chr2pos2sv_info;
map<string, map<int,  sv_var_info> > virus2pos2sv_info;
map<string, map<int,  insertion_info> > chr2pos2insertion_info;
map<int,  var_info> pos2var_info;
map<int, sv_var_info> pos2sv_info;
map<string, map<int, int> > chr2nbg2nend;
map<pair<int, int>, pair<int,char> > pos2allele;
map<pair<int, int>, map<pair<int, int>, pair<int,char> > > pos2cppos_allele;
map<string, string> virus2seq;
map<string, string> chr2seq;
vector<string> vec_chr;
vector<string> vec_virus;
pair<map<string, string>::iterator, bool> ret;
//test code
int total_len1 = 0, total_len2 = 0;
//test code

MVSCOption mvscoption;
variation_ratio var;
double ins_len_ratio[20];
double del_len_ratio[20];
ofstream FOUT_B, FOUTSNP_B, FOUTINDEL_B, FOUTSV_B, FOUTSEQ_B, FOUTCNV_B, FOUT_P, FOUTSNP_P, FOUTINDEL_P, FOUTSV_P, FOUTSEQ_P, FOUTCNV_P, FOUT_C, FOUTSNP_C, FOUTINDEL_C, FOUTSV_C, FOUTSEQ_C, FOUTCNV_C;

int OutputChrWithoutVariation(string chr, string hg_seq);
int Simulate(string, string);
void get_ins_ratio()
{
		double a=0.53,b=1.60;
		int i;
		for(i=0;i<20;i++)
		{
				if(i>0)
				{
						ins_len_ratio[i]=ins_len_ratio[i-1]+a*pow(i+1,-1*b);
				}
				else
				{
						ins_len_ratio[i]=a*pow(i+1,-1*b);
				}
		}
}
void get_del_ratio()
{
		double a=0.48,b=1.51;
		int i;
		for(i=0;i<20;i++)
		{
				if(i>0)
				{
						del_len_ratio[i]=del_len_ratio[i-1]+a*pow(i+1,-1*b);
				}
				else
				{
						del_len_ratio[i]=a*pow(i+1,-1*b);
				}
		}
}
int main(int argc, char *argv[])
{
	int opt=0;
	static const char*optString="o:v:d:D:i:l:s:S:V:N:n:b:p:t:c:?"; //option parameters
	opt=getopt(argc,argv,optString);
	while(opt!=-1)
	{
		switch(opt)
		{
		case 'o':
			mvscoption.output_dir=optarg;
			break;
		case 'v':
			mvscoption.simu_str=optarg;
			break;
		case 'd':
			mvscoption.hg_ref=optarg;
			break;
//		case 'h':
//			mvscoption.chr=optarg;
//			break;
		case 'D':
			mvscoption.virus_ref=optarg;
			break;
		case 's':
		    mvscoption.snp_variation_ratio=atof(optarg);
			break;
		case 'i':
		    mvscoption.indel_variation_ratio=atof(optarg);
			break;
		case 'S':
		    mvscoption.insertion_rate=atof(optarg);
			break;
		case 'V':
		    mvscoption.inversion_rate=atof(optarg);
			break;
		case 'N':
		    mvscoption.deletion_rate=atof(optarg);
			break;
		case 'l':
		    mvscoption.sv_variation_ratio=atof(optarg);
			break;
		case 't':
		    mvscoption.cancer_sv_variation_ratio=atof(optarg);
			break;
		case 'n':
			mvscoption.narea_file = optarg;
			break;
		case 'b':
			mvscoption.blood_variation_ratio=atof(optarg);
			break;
		case 'p':
			mvscoption.para_variation_ratio=atof(optarg);
			break;
		case 'c':
			mvscoption.cancer_variation_ratio=atof(optarg);
			break;
		case'?':
			display_usage();
			break;
		default:
			break;
		}
		opt=getopt(argc,argv,optString);
	}
	if(mvscoption.output_dir.empty() || mvscoption.hg_ref.empty())
	{
		display_usage();
		exit(EXIT_FAILURE);
	}
	for(int i=0;i<5;i++)
	{
		switch(mvscoption.simu_str[i])
		{
			case '1':
				break;
			case '0':
				VAR_type[i]=false;
				break;
		}
	}
	get_del_ratio();
	get_ins_ratio();
	
	mvscoption.Compute();
	OpenFiles(mvscoption, FOUT_B, FOUTSNP_B, FOUTINDEL_B, FOUTSV_B, FOUTSEQ_B, FOUTCNV_B, FOUT_P, FOUTSNP_P, FOUTINDEL_P, FOUTSV_P, FOUTSEQ_P, FOUTCNV_P, FOUT_C, FOUTSNP_C, FOUTINDEL_C, FOUTSV_C, FOUTSEQ_C, FOUTCNV_C);

	ifstream fin;
	string virus_seq("");//first sequence
	if(VAR_type[2] && !mvscoption.virus_ref.empty())//long insertion
	{
		ReadVirusSeq(mvscoption.virus_ref, virus2seq, vec_virus);
	}
	else if(VAR_type[2])
	{
		if(mvscoption.virus_ref.empty()) {cout <<"didn't input virus sequence file, do not simulate exogenous insertion"<<endl;}
	}

	if (VAR_type[2] || VAR_type[3] || VAR_type[4])
	{
		if (mvscoption.narea_file.empty())
		{
			cout << "Warning:didn't input the file which record N area, the SV and CNV may happen in N area" << endl;
		}
		else 
		{
			fin.open(mvscoption.narea_file.c_str());
			if (!fin) {cerr << "Cannot open the file which record N area." << endl;exit(EXIT_FAILURE);}
			string chrN, temp;
			int start = 0, end = 0;
			map<string, map<int, int> >::iterator chr2nbg2nend_it;	
			while (fin >> chrN >> start >> end)
			{
				getline(fin, temp);
				chr2nbg2nend_it = chr2nbg2nend.find(chrN);
				if (chr2nbg2nend_it == chr2nbg2nend.end()) {
					map<int, int> nbg2nend;
					nbg2nend.insert(make_pair(start, end));
					chr2nbg2nend.insert(make_pair(chrN, nbg2nend));
				}
				else {
					chr2nbg2nend_it->second.insert(make_pair(start, end));
				}
			}
			fin.close();
		}
	}
	//variation ratio Initialization
	var.compute(mvscoption);


	fin.open(mvscoption.hg_ref.c_str());
	if(!fin.is_open())
	{
		cerr<<"Please check whether this file exists or the file name is misspelled!"<<endl;
		exit(EXIT_FAILURE);
	}

	ifstream FIN(mvscoption.hg_ref.c_str());	
	if (!FIN) { cerr << "Open file " << mvscoption.hg_ref << " failure!" << endl; exit(1); }
	string chr_or_seq;
	string chr, hg_seq;
	while (getline(FIN, chr_or_seq)) {
		if (chr_or_seq[0] == '>') {
			if (!chr.empty()) {
				if (hg_seq.size()<=1000) {
					OutputChrWithoutVariation(chr, hg_seq);
				}
				else {
					transform(hg_seq.begin(), hg_seq.end(), hg_seq.begin(), ::toupper);
					ret = chr2seq.insert(make_pair(chr, hg_seq));
					if (ret.second) {
						vec_chr.push_back(chr);
					}
					else {
						cerr << "Warning: " << chr << " repeated in the reference file " << mvscoption.hg_ref << endl;
					}
				}
				hg_seq.clear();
			}
			chr = chr_or_seq.substr(1);
			chr = GetStringBeforeSpace(chr);
		}
		else {
			hg_seq.append(chr_or_seq);
		}
	}
	if (hg_seq.size()<=1000) {
		OutputChrWithoutVariation(chr, hg_seq);
	}
	else {
		transform(hg_seq.begin(), hg_seq.end(), hg_seq.begin(), ::toupper);
		ret = chr2seq.insert(make_pair(chr, hg_seq));
		if (ret.second) {
			vec_chr.push_back(chr);
		}
		else {
			cerr << "Warning: " << chr << " repeated in the reference file " << mvscoption.hg_ref << endl;
		}
	}
	hg_seq.clear();


	srand((unsigned)time(NULL)); //Blood
	//Simulation_start(chr2seq, vec_chr, var, virus2seq, chr2pos2var_info, chr2pos2sv_info, chr2nbg2nend, pos2allele, pos2cppos_allele,'B');
	Simulation_start(chr2seq, vec_chr, virus2seq, vec_virus, chr2pos2var_info, chr2pos2sv_info, virus2pos2sv_info, chr2pos2insertion_info, chr2nbg2nend, var, 'B');
	cerr << "Germline:Variation simulation finished!" << endl;
	VariationIntegration(FOUT_B, FOUTSNP_B, FOUTINDEL_B, FOUTSV_B, chr2seq, virus2seq, vec_chr, vec_virus, chr2pos2var_info, chr2pos2insertion_info, chr2pos2sv_info, virus2pos2sv_info, chr2nbg2nend, 'B');
	Simulation_start(chr2seq, vec_chr, virus2seq, vec_virus, chr2pos2var_info, chr2pos2sv_info, virus2pos2sv_info, chr2pos2insertion_info, chr2nbg2nend, var, 'P');
	cerr << "Para:Variation simulation finished!" << endl;
	VariationIntegration(FOUT_P, FOUTSNP_P, FOUTINDEL_P, FOUTSV_P, chr2seq, virus2seq, vec_chr, vec_virus, chr2pos2var_info, chr2pos2insertion_info, chr2pos2sv_info, virus2pos2sv_info, chr2nbg2nend, 'B');
	Simulation_start(chr2seq, vec_chr, virus2seq, vec_virus, chr2pos2var_info, chr2pos2sv_info, virus2pos2sv_info, chr2pos2insertion_info, chr2nbg2nend, var, 'C');
	cerr << "Cancer:Variation simulation finished!" << endl;
	VariationIntegration(FOUT_C, FOUTSNP_C, FOUTINDEL_C, FOUTSV_C, chr2seq, virus2seq, vec_chr, vec_virus, chr2pos2var_info, chr2pos2insertion_info, chr2pos2sv_info, virus2pos2sv_info, chr2nbg2nend, 'B');

//	Simulation_start(chr2seq, vec_chr, virus2seq, vec_virus, chr2pos2var_info, chr2pos2sv_info, virus2pos2sv_info, chr2pos2insertion_info, chr2nbg2nend, var, 'P');
//	cerr << "Para:Variation simulation finished!" << endl;
//	VariationIntegration(FOUT_P, chr2seq, virus2seq, vec_chr, vec_virus, chr2pos2var_info, chr2pos2insertion_info, chr2pos2sv_info, virus2pos2sv_info, chr2nbg2nend, 'P');
//
//	Simulation_start(chr2seq, vec_chr, virus2seq, vec_virus, chr2pos2var_info, chr2pos2sv_info, virus2pos2sv_info, chr2pos2insertion_info, chr2nbg2nend, var, 'C');
//	cerr << "Cancer:Variation simulation finished!" << endl;
//	VariationIntegration(FOUT_C, chr2seq, virus2seq, vec_chr, vec_virus, chr2pos2var_info, chr2pos2insertion_info, chr2pos2sv_info, virus2pos2sv_info, chr2nbg2nend, 'C');

	//	variation_integration(hg_seq, chr2pos2var_info_it->second, chr2pos2sv_info_it->second, outseq1, outseq2);
		//Output(FOUT_P, FOUTSNP_P, FOUTINDEL_P, FOUTSV_P, FOUTSEQ_P, FOUTCNV_P, chr, outseq1, outseq2, chr2pos2var_info_it->second, chr2pos2sv_info_it->second, pos2allele, pos2cppos_allele);

	return 0;
}

int OutputChrWithoutVariation(string chr, string hg_seq)
{
	cerr << "The chr's length must be longer than 1000, do not simulate mutation in chromosome " << chr << endl;
	FOUT_B << ">" << chr << "/1" << endl;
	OutputSeq(FOUT_B, hg_seq);
	FOUT_B << ">" << chr << "/2" << endl;
	OutputSeq(FOUT_B, hg_seq);
	FOUT_P << ">" << chr << "/1" << endl;
	OutputSeq(FOUT_P, hg_seq);
	FOUT_P << ">" << chr << "/2" << endl;
	OutputSeq(FOUT_P, hg_seq);
	FOUT_C << ">" << chr << "/1" << endl;
	OutputSeq(FOUT_C, hg_seq);
	FOUT_C << ">" << chr << "/2" << endl;
	OutputSeq(FOUT_C, hg_seq);
}

int Simulate(string chr, string hg_seq)
{
//	srand((unsigned)time(NULL)); //Blood
//	string outseq1, outseq2;
//	Simulation_start(chr, hg_seq, var, virus2seq, chr2pos2var_info,pos2sv_info,nbg2nend,pos2allele,pos2cppos_allele,'B');
//		variation_integration(hg_seq,pos2var_info,pos2sv_info,outseq1,outseq2);
	
		//Output(FOUT_B, FOUTSNP_B, FOUTINDEL_B, FOUTSV_B, FOUTSEQ_B, FOUTCNV_B, chr, outseq1, outseq2, pos2var_info, pos2sv_info, pos2allele, pos2cppos_allele);
		//para
//		outseq1 = "";
//		outseq2 = "";
//		Simulation_start(chr, hg_seq,var, virus2seq, chr2pos2var_info,pos2sv_info,nbg2nend,pos2allele,pos2cppos_allele,'P');
//		variation_integration(hg_seq,pos2var_info,pos2sv_info,outseq1,outseq2);
//		Output(FOUT_P, FOUTSNP_P, FOUTINDEL_P, FOUTSV_P, FOUTSEQ_P, FOUTCNV_P, chr, outseq1, outseq2, pos2var_info, pos2sv_info, pos2allele, pos2cppos_allele);
		//Cancer
//		outseq1 = "";
//		outseq2 = "";
//		Record_deletion(mvscoption.retain_paraTocancer_ratio,pos2var_info,pos2sv_info,pos2allele,pos2cppos_allele);//retain part of para variation
//		Simulation_start(chr, hg_seq, var, virus2seq, chr2pos2var_info,pos2sv_info,nbg2nend,pos2allele,pos2cppos_allele,'C');
//		variation_integration(hg_seq,pos2var_info,pos2sv_info,outseq1,outseq2);
//		Output(FOUT_C, FOUTSNP_C, FOUTINDEL_C, FOUTSV_C, FOUTSEQ_C, FOUTCNV_C, chr, outseq1, outseq2, pos2var_info, pos2sv_info, pos2allele, pos2cppos_allele);

//		pos2var_info.clear();
//		pos2sv_info.clear();
//		pos2allele.clear();
//		pos2cppos_allele.clear();
	return 0;
}
