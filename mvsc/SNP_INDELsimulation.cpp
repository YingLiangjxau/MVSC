#include "cgvsimulation.h"
#include <vector>

#define VERSION "2.0"
#define REVISION 6 
void snp_simulation(const string &chr, const string& seq, double snp_variation_ratio, map<string, map<int, var_info> >& chr2pos2var_info, char level)
{
	map<string, map<int, var_info> >::iterator chr2pos2var_info_it = chr2pos2var_info.find(chr);
	if (chr2pos2var_info_it == chr2pos2var_info.end()) {
		map<int, var_info> pos2var_info;
		chr2pos2var_info_it = chr2pos2var_info.insert(make_pair(chr, pos2var_info)).first;
	}

	int Chr_len=seq.length();//length=133,851,895  10*8
	int sum=int(Chr_len*snp_variation_ratio);
	cout<<"Initialization SNP count="<<sum<<endl;
	int times=0, ref_position=-1;
	while (sum)
	{
	   	//int count=0;
	    ref_position=rand()%Chr_len;
		//  cerr<<"SNP count="<<++count<<endl;
		if((seq.at(ref_position)!='N')&&(ref_position-flank>=0)&&(ref_position+flank<=Chr_len))
		{
			if (isCanInsert(ref_position + 1, 1, chr2pos2var_info_it->second))
			{
				var_info ins;
			 	//snp sites
				ins.ref = seq[ref_position];//string start from 0
				string AB=SNP_var(ins.ref);
				//storage structure
				//var_info ins = {"snp", AB[0], AB[1], 1, 1, level};
				ins.var_type = "snp";
				ins.seq1 = AB[0];
				ins.length1 = 1;
				ins.seq2 = AB[1];
				ins.length2 = 1;
				ins.period = level;
				//var_info ins = {"snp", AB[0], AB[1], 1, 1, level};
				chr2pos2var_info_it->second.insert(make_pair(ref_position+1, ins));//map<int, var_info>
				++times;
			 	if(times%1000==0&&times!=0)
			 	{
					cout<<times<<" SNP simulation finished ..."<<endl;
				}
				if(times==sum)
				{
					cout<<times<<" SNP simulation finished!"<<endl;
					break;
				}

			}
		}
	}
}

//small indel seq length is less than 20, and 90 percent is less than 5, without N area 
//assume homozygous vs heterozygous 1:4 for small indel
void indel_simulation(const string &chr, const string& seq, double indel_variation_ratio, map<string, map<int, var_info> >& chr2pos2var_info, char level)
{	
	map<string, map<int, var_info> >::iterator chr2pos2var_info_it = chr2pos2var_info.find(chr);
	if (chr2pos2var_info_it == chr2pos2var_info.end()) {
		map<int, var_info> pos2var_info;
		chr2pos2var_info_it = chr2pos2var_info.insert(make_pair(chr, pos2var_info)).first;
	}

	int Chr_len=seq.length();
	int sum=int(Chr_len*indel_variation_ratio);	//if sum=10000 times
	int del_count=0;
	int ins_count=0;
	int del_Homozygous_count=0;
	int ins_Homozygous_count=0;
	int del_Heterozygous_count=0;
	int ins_Heterozygous_count=0;
	int del_Heterozygous_chainA=0;
	int ins_Heterozygous_chainA=0;
	int del_Heterozygous_chainB=0;
	int ins_Heterozygous_chainB=0;


	//variation type
	for(int i=0;i<sum;i++)
	{
		int insOrdel=rand()%2;
		if(insOrdel==0){del_count++;}
	}
	ins_count=sum-del_count;
	//#pragma omp critical
	//cout<<"Initialization insertion:\tins_count="<<ins_count<<"\tdel_count="<<del_count<<endl;
	//deletion Homozygous OR Heterozygous
	for(int i=0;i<del_count;i++)
	{
		int zygous=rand()%5;//homozygous:heterozygous=1:4
		if(zygous==0){del_Homozygous_count++;}
	}
	//insertion Homozygous OR Heterozygous
	for(int i=0;i<ins_count;i++)
	{
		int zygous=rand()%5;//homozygous:heterozygous=1:4
		if(zygous==0){ins_Homozygous_count++;}
	}
	del_Heterozygous_count=del_count-del_Homozygous_count;
    ins_Heterozygous_count=ins_count-ins_Homozygous_count;

/*#pragma omp critical
{
	cout<<"Initialization deletion:\tdel_Heterozygous_count="<<del_Heterozygous_count<<"\tdel_Homozygous_count="<<del_Homozygous_count<<endl;
	cout<<"Initialization insertion:\tins_Heterozygous_count="<<ins_Heterozygous_count<<"\tins_Homozygous_count="<<ins_Homozygous_count<<endl;
}*/
	//deletion		ChainAB

	for(int i=0;i<del_Heterozygous_count;i++)
	{
		int zygous=rand()%2;
		if(zygous==0){del_Heterozygous_chainA++;}
	}
	//insertion	ChainAB
	for(int i=0;i<ins_Heterozygous_count;i++)
	{
		int zygous=rand()%2;
		if(zygous==0){ins_Heterozygous_chainA++;}
	}
	del_Heterozygous_chainB=del_Heterozygous_count-del_Heterozygous_chainA;
	ins_Heterozygous_chainB=ins_Heterozygous_count-ins_Heterozygous_chainA;
//#pragma omp critical
//{
	cout<<"Initialization insertion:\tins_count="<<ins_count<<"\tdel_count="<<del_count<<endl;
	cout<<"Initialization deletion:\tdel_Heterozygous_count="<<del_Heterozygous_count<<"\tdel_Homozygous_count="<<del_Homozygous_count<<endl;
	cout<<"Initialization insertion:\tins_Heterozygous_count="<<ins_Heterozygous_count<<"\tins_Homozygous_count="<<ins_Homozygous_count<<endl;
	cout<<"Initialization deletion:\tdel_Heterozygous_chainB="<<del_Heterozygous_chainB<<"\tdel_Heterozygous_chainA="<<del_Heterozygous_chainA<<endl;
	cout<<"Initialization insertion:\tins_Heterozygous_chainB="<<ins_Heterozygous_chainB<<"\tins_Heterozygous_chainA="<<ins_Heterozygous_chainA<<endl;
//}

	int times=0;
	while(del_count)//deletion Variation times
	{
		int lens=get_del_length();
		//	cout<<"lens="<<lens<<endl;
		int ref_position=-1;
		int found=1; 
		string str("");
		//int count=0;

		ref_position=rand()%Chr_len;
	    if(ref_position+lens+flank>Chr_len || ref_position-flank<0){continue;}//deletion is out of Reference

		//cerr<<"deletion count="<<++count<<endl;
		str=seq.substr(ref_position,lens);
		found=str.find('N',0);
		if(found==string::npos)
		{
			if (!isCanInsert(ref_position + 1, lens, chr2pos2var_info_it->second)) {continue;}
			//isSmallAndSmallOverlap(ref_position + 1, lens, low_bound_map_it, up_bound_map_it, pos2var_info.end())){continue;}
			//Homozygous and Heterozygous for Chain AB
			var_info ins;
			ins.period=level;
			ins.var_type = "del";
			if(times<=del_Homozygous_count)	//Homozygous
			{	
				ins.seq1=str;
				ins.length1=lens;
				ins.seq2=ins.seq1;
				ins.length2=lens;
			}//Heterozygous chainA
			else if(times<=del_Homozygous_count+del_Heterozygous_chainA)
			{				
				ins.seq1=str;
				ins.length1=lens;
				ins.seq2="H";
				ins.length2=0;
			}//Heterozygous chainB
			else if(times<=del_count)
			{
				ins.seq1="H";
				ins.length1=0;						
				ins.seq2=str;
				ins.length2=lens;
			}
			chr2pos2var_info_it->second.insert(make_pair(ref_position+1,ins));//map<int, pair<string, var_info> >
			++times;

			if(times%1000==0&&times!=0)
			{
			 	//#pragma omp critical
				cout<<times<<" small Deletion Simulation finished ..."<<endl;
			}
			if(times == del_count)
			{
			 	//#pragma omp critical
				cout << times << " small Deletion Simulation finished!" << endl;
				break;
			}
		}

	}

	times=0;
	while (ins_count)
	{
		int lens=get_ins_length();
		//int count=0;
		int ref_position=rand()%Chr_len;
		if (ref_position-flank<0 || ref_position+1+flank>Chr_len){continue;}//insertion is out of Reference
		//	cerr<<"deletion count="<<++count<<endl;
	 	if(seq.at(ref_position)!='N')
		{ 
			if (!isCanInsert(ref_position + 1, 1, chr2pos2var_info_it->second)) {continue;}
			//Select the insertion sequence
			string fragment("");//insertion sequence
			{
				const char CCH[] = "ATCG";
				char ch[21]={0};
				for (int i = 0; i < lens; ++i)
				{
					int x = rand() / (RAND_MAX / (sizeof(CCH) - 1));
					ch[i] = CCH[x];
				}
				fragment.assign(ch);
			}
			var_info ins;
			if(times<=ins_Homozygous_count)//Homozygous
			{	
				ins.seq1=fragment;
				ins.length1=lens;
				ins.seq2=fragment;
				ins.length2=lens;
			}
			//Heterozygous chainA
			else if(times<=ins_Homozygous_count+ins_Heterozygous_chainA)
			{ 
				ins.seq1=fragment;
				ins.length1=lens;
				ins.seq2="H";
				ins.length2=0;
			}
			//Heterozygous chainB
			else if(times<=ins_count)
			{									 
				ins.seq1="H";
				ins.length1=0;						
				ins.seq2=fragment;
				ins.length2=lens;
			}
			ins.period=level;
			ins.var_type="ins";
			chr2pos2var_info_it->second.insert(make_pair(ref_position+1,ins));
			++times;
			
			//test code
			cout<<times<<" Test:small Insertion Simulation ..."<<endl;
			//test code

			if(times%1000==0&&times!=0)
			{
				//#pragma omp critical
				cout<<times<<" small Insertion Simulation finished ..."<<endl;
			}
			if(times==ins_count)
			{
				cout<<times<<" small Insertion Simulation finished ..."<<endl;
				break;
			}
		}
	}

}


