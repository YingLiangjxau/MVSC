MVSC (Multi-Variation Simulator of Cancer genome), a tool to simulate common genome variants involving single nucleotide polymorphisms, small insertion and deletion polymorphisms, and structural variations, which are analogous to human somatically acquired variations. Three sets of variations embedded genomic sequences in different periods are dynamically and sequentially simulated one by one. In cancer genome simulation, complex SV is imported because this variation helps elucidate tumor genome structure. Variations with different sizes can coexist in the genome region, which are not mutually exclusive.

The software is now compatible with linux operation system. The example folder is a toy example,and source code is in mvsc folder, you can directly use the mvsc (exe file) to simulate variations or compile source code according to makefile.

MVSC:&ensp;a multi-variation simulator of cancer genome
<br>Version 2.0
<br>Revision 6
<br>Usage:	mvsc
<br>Required parameters:
<br>&emsp;&emsp;	-o	<str>	output directory
<br>&emsp;&emsp;	-d	<str>	input human reference

<br>Optional:
<br>&emsp;&emsp;	-v	<str>	five number string, it  represent in order[snp, small indels, large insersion, large deletion, inversion], 0 means Not simulate this type of variation,1 means simulation(default:11111)
<br>&emsp;&emsp;	-D	<str>	input exogenous sequence for long indel(default:all of the large insertion are endogenous)
<br>&emsp;&emsp;	-b	<double>	the variation ratio in germline(default:0.0006)
<br>&emsp;&emsp;	-s	<double>	the snp ratio based on variation number in germline(default:0.877)
<br>&emsp;&emsp;	-i	<double>	the indel ratio based on variation number in germline(default:0.105)
<br>&emsp;&emsp;	-l	<double>	the sv ratio based on variation number in germline(default:0.002)
<br>&emsp;&emsp;	-p	<double>	the variation ratio in para based on germline(default:0.00001)
<br>&emsp;&emsp;	-S	<double>	the insertion ratio among SV types(default:0.3)
<br>&emsp;&emsp;	-V	<double>	the inversion ratio among SV types(default:0.25)
<br>&emsp;&emsp;	-N	<double>	the deletion ratio among SV types(default:0.45)
<br>&emsp;&emsp;	-c	<double>	the variation ratio in cancer based on para(default:0.00003)
<br>&emsp;&emsp;	-t	<double>	the sv ratio based on variation number in cancer,which should be larger than the parametre set -l(default:10*l)

<br>If you have any questions,please tell us: aliang1229@126.com
<br>The MVSC software package on https://github.com/qiukunlong/variation_simulation is the earliest version. And the newest version of MVSC was on https://github.com/YingLianghnu/MVSC.

