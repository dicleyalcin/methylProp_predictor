from Bio import motifs 
from joblib import Parallel, delayed 
import multiprocessing 
import os 

num_cores= multiprocessing.cpu_count() 

# Dataset variables
dataset = "Feltus"
TFBS_database = "UniPROBE" 
dataset_property = "unbalanced" 

inputs = range(500)

def processInput(i):     

	#------------------------------------------------------------------------------------------------------------------
	# # # # # # #  ------- Create split directories and put all .fa files in each split ----------- # #------------------------------------------------------------------------------------------------------------------
	
	#return os.system("mkdir /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d && mkdir /home/data/dicleyalcin/%s_%s/MEME/TrainR/TrainR_Split_%d && mkdir /home/data/dicleyalcin/%s_%s/MEME/TestP/TestP_Split_%d && mkdir /home/data/dicleyalcin/%s_%s/MEME/TestR/TestR_Split_%d && mv /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d.fa /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/ && mv /home/data/dicleyalcin/%s_%s/MEME/TrainR/TrainR_Split_%d.fa /home/data/dicleyalcin/%s_%s/MEME/TrainR/TrainR_Split_%d/ && mv /home/data/dicleyalcin/%s_%s/MEME/TestP/TestP_Split_%d.fa /home/data/dicleyalcin/%s_%s/MEME/TestP/TestP_Split_%d/ && mv /home/data/dicleyalcin/%s_%s/MEME/TestR/TestR_Split_%d.fa /home/data/dicleyalcin/%s_%s/MEME/TestR/TestR_Split_%d/" % (dataset,dataset_property,(i+1),dataset,dataset_property,(i+1),dataset,dataset_property,(i+1),dataset,dataset_property,(i+1),dataset,dataset_property,(i+1),dataset,dataset_property,(i+1),dataset,dataset_property,(i+1),dataset,dataset_property,(i+1),dataset,dataset_property,(i+1),dataset,dataset_property,(i+1),dataset,dataset_property,(i+1),dataset,dataset_property,(i+1)))

	#------------------------------------------------------------------------------------------------------------------
	# # # # # # # # # # # ----- Running psp-gen on all Training Splits  ----------- # # # # # # # # # # # 
	#------------------------------------------------------------------------------------------------------------------

	#return os.system("cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d && psp-gen -pos TrainP_Split_%d.fa -neg ../../TrainR/TrainR_Split_%d/TrainR_Split_%d.fa > primary_prones_split_%d.psp && psp-gen -pos ../../TrainR/TrainR_Split_%d/TrainR_Split_%d.fa -neg TrainP_Split_%d.fa > primary_resistants_split_%d.psp" % (dataset,dataset_property,(i+1),(i+1),(i+1),(i+1),(i+1),(i+1),(i+1),(i+1),(i+1))) 


	#------------------------------------------------------------------------------------------------------------------
	# # # # # # # # # # # ----- Running MEME on all Training Splits (FIND PMs for 500 splits, and RMs for 500 splits)  
	#------------------------------------------------------------------------------------------------------------------


	#return os.system("cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d && meme TrainP_Split_%d.fa -psp primary_prones_split_%d.psp  -dna -revcomp -oc meme_out_prones_split_%d -mod zoops -nmotifs 30 -minw 6 -maxw 50 -maxsize 500000 && meme ../../TrainR/TrainR_Split_%d/TrainR_Split_%d.fa -psp primary_resistants_split_%d.psp -dna -revcomp -oc meme_out_resistants_split_%d -mod zoops -nmotifs 30 -minw 6 -maxw 50 -maxsize 500000" % (dataset,dataset_property,(i+1),(i+1),(i+1),(i+1),(i+1),(i+1),(i+1),(i+1)))																								



	#------------------------------------------------------------------------------------------------------------------
	# # # # # # # # # # # # # # # Extracting consensus sequences of motifs found (for each split) ----- # # # # # # # 
	#------------------------------------------------------------------------------------------------------------------

	#return os.system("PATH=$HOME/meme/bin:$PATH &&  cd MEME/TrainP/TrainP_Split_%d && cd MEME_outputs && python /home/dyalcin/Scripts/extractor.py > consensus_motifs_split_%d.txt &&  cd ../../../../MEME/TrainR/TrainR_Split_%d && cd MEME_outputs && python /home/dyalcin/Scripts/extractor.py > consensus_motifs_split_%d.txt" % ((i+1),(i+1),(i+1),(i+1)))


	#------------------------------------------------------------------------------------------------------------------
 	# # # # # # # # # # # # # #Extracting XML files of MEME from each split to a separate folder # # # # # # # # # #
 	#------------------------------------------------------------------------------------------------------------------


	#return os.system("PATH=$HOME/meme/bin:$PATH && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_prones_split_%d/ && mv meme.xml meme_trainP_split_%d.xml && cp meme_trainP_split_%d.xml /home/data/dicleyalcin/%s_%s/MEME_significant_motifs/ && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_resistants_split_%d/ && 	mv meme.xml meme_trainR_split_%d.xml && cp meme_trainR_split_%d.xml /home/data/dicleyalcin/%s_%s/MEME_significant_motifs/" % (dataset,dataset_property,(i+1),(i+1),(i+1),(i+1),dataset,dataset_property,dataset,dataset_property,(i+1),(i+1),(i+1),(i+1),dataset,dataset_property))



	#------------------------------------------------------------------------------------------------------------------
	# # # # # # # # # #  Extracting significant motifs in each spliit using the xml files of MEMEs # # # # # # # #     
	#------------------------------------------------------------------------------------------------------------------


	#return os.system("PATH=$HOME/meme/bin:$PATH && cd /home/data/dicleyalcin/%s_%s/MEME_significant_motifs && python /home/dyalcin/Scripts/meme_sigMotif_extractor.py" % ((dataset, dataset_property)))    


	#------------------------------------------------------------------------------------------------------------------
	# # # # # # # #------ MAST runs for TRAIN and TEST data using MEME motifs created for each split 
	#------------------------------------------------------------------------------------------------------------------

	"""
	return os.system("PATH=$HOME/meme/bin:$PATH && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_prones_split_%d && mast meme.txt -oc MAST_outputs_MotifProne_dbProne ../TrainP_Split_%d.fa -remcorr -ev 100000 && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_prones_split_%d && mast meme.txt -oc MAST_outputs_MotifProne_dbResistant ../../../TrainR/TrainR_Split_%d/TrainR_Split_%d.fa -remcorr -ev 100000 && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_resistants_split_%d && mast meme.txt -oc MAST_outputs_MotifResistant_dbProne ../TrainP_Split_%d.fa -remcorr -ev 100000 && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_resistants_split_%d && mast meme.txt -oc MAST_outputs_MotifResistant_dbResistant ../../../TrainR/TrainR_Split_%d/TrainR_Split_%d.fa -remcorr -ev 100000 && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_prones_split_%d && mast meme.txt -oc MAST_outputs_MotifProne_dbProne_TEST ../../../TestP/TestP_Split_%d/TestP_Split_%d.fa -remcorr -ev 100000 && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_prones_split_%d && mast meme.txt -oc MAST_outputs_MotifProne_dbResistant_TEST ../../../TestR/TestR_Split_%d/TestR_Split_%d.fa -remcorr -ev 100000 && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_resistants_split_%d && mast meme.txt -oc MAST_outputs_MotifResistant_dbProne_TEST ../../../TestP/TestP_Split_%d/TestP_Split_%d.fa -remcorr -ev 100000 && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_resistants_split_%d && mast meme.txt -oc 					MAST_outputs_MotifResistant_dbResistant_TEST ../../../TestR/TestR_Split_%d/TestR_Split_%d.fa -remcorr -ev 100000" % (dataset,dataset_property, (i+1),(i+1),(i+1),dataset, dataset_property, (i+1),(i+1),(i+1),(i+1), dataset, dataset_property,(i+1),(i+1),(i+1),dataset,dataset_property,(i+1),(i+1),(i+1),(i+1),dataset,dataset_property,(i+1),(i+1),(i+1),(i+1),dataset,dataset_property,(i+1),(i+1),(i+1),(i+1),dataset,dataset_property,(i+1),(i+1),(i+1),(i+1),dataset,dataset_property,(i+1),(i+1),(i+1),(i+1)))					    
	
	"""
	#------------------------------------------------------------------------------------------------------------------
	# # # # # # # #------ MAST runs for TRAIN and TEST data using MEME motifs created for each split (overfit case)
	#------------------------------------------------------------------------------------------------------------------

	
	#return os.system("PATH=$HOME/meme/bin:$PATH && 																									cd /home/dyalcin/data/%s_%s/MEME-overfit/ && mast meme_prones.xml -oc MAST_outputs_MotifProne_dbProne_split_%d ../MEME/TrainP/TrainP_Split_%d/TrainP_Split_%d.fa -remcorr -ev 100000 && cd /home/dyalcin/data/%s_%s/MEME-overfit/ && mast meme_prones.xml -oc MAST_outputs_MotifProne_dbResistant_split_%d ../MEME/TrainR/TrainR_Split_%d/TrainR_Split_%d.fa -remcorr -ev 100000 && cd /home/dyalcin/data/%s_%s/MEME-overfit/ && mast meme_resistants.xml -oc MAST_outputs_MotifResistant_dbProne_split_%d ../MEME/TrainP/TrainP_Split_%d/TrainP_Split_%d.fa -remcorr -ev 100000 && cd /home/dyalcin/data/%s_%s/MEME-overfit/ && mast meme_resistants.xml -oc MAST_outputs_MotifResistant_dbResistant_split_%d ../MEME/TrainR/TrainR_Split_%d/TrainR_Split_%d.fa -remcorr -ev 100000 && cd /home/dyalcin/data/%s_%s/MEME-overfit/ && mast meme_prones.xml -oc MAST_outputs_MotifProne_dbProne_TEST_split_%d ../MEME/TestP/TestP_Split_%d/TestP_Split_%d.fa -remcorr -ev 100000 && cd /home/dyalcin/data/%s_%s/MEME-overfit/ && mast meme_prones.xml -oc MAST_outputs_MotifProne_dbResistant_TEST_split_%d ../MEME/TestR/TestR_Split_%d/TestR_Split_%d.fa -remcorr -ev 100000 && cd /home/dyalcin/data/%s_%s/MEME-overfit/ && mast meme_resistants.xml -oc MAST_outputs_MotifResistant_dbProne_TEST_split_%d ../MEME/TestP/TestP_Split_%d/TestP_Split_%d.fa -remcorr -ev 100000 && cd /home/dyalcin/data/%s_%s/MEME-overfit/ && mast meme_resistants.xml -oc MAST_outputs_MotifResistant_dbResistant_TEST_split_%d ../MEME/TestR/TestR_Split_%d/TestR_Split_%d.fa -remcorr -ev 100000" % (dataset, dataset_property,  (i+1),(i+1),(i+1),dataset,dataset_property, (i+1), (i+1),(i+1), dataset,dataset_property, (i+1),(i+1),(i+1), dataset,dataset_property, (i+1),(i+1),(i+1), dataset,dataset_property, (i+1),(i+1),(i+1), dataset,dataset_property, (i+1),(i+1),(i+1), dataset,dataset_property, (i+1),(i+1),(i+1), dataset,dataset_property, (i+1),(i+1),(i+1)))
	

	#------------------------------------------------------------------------------------------------------------------
	# # # # # # # # # # # # # # #------ unique MAST xml files are created and copied in MAST  ----- # # # # # # # # # 
	#------------------------------------------------------------------------------------------------------------------
	"""
	return os.system("PATH=$HOME/meme/bin:$PATH && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_prones_split_%d/MAST_outputs_MotifProne_dbProne &&  mv mast.xml mast_PP_Split_%d.xml &&  cp mast_PP_Split_%d.xml ../../../../../MAST && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_prones_split_%d/MAST_outputs_MotifProne_dbResistant && mv mast.xml mast_PR_Split_%d.xml &&  cp mast_PR_Split_%d.xml ../../../../../MAST && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_resistants_split_%d/MAST_outputs_MotifResistant_dbProne &&  mv mast.xml mast_RP_Split_%d.xml &&  cp mast_RP_Split_%d.xml ../../../../../MAST && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_resistants_split_%d/MAST_outputs_MotifResistant_dbResistant &&  mv mast.xml mast_RR_Split_%d.xml && cp mast_RR_Split_%d.xml ../../../../../MAST && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_prones_split_%d/MAST_outputs_MotifProne_dbProne_TEST && mv mast.xml mast_PP_Split_%d_TEST.xml && cp mast_PP_Split_%d_TEST.xml ../../../../../MAST && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_prones_split_%d/MAST_outputs_MotifProne_dbResistant_TEST &&  mv mast.xml mast_PR_Split_%d_TEST.xml && cp mast_PR_Split_%d_TEST.xml ../../../../../MAST && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_resistants_split_%d/MAST_outputs_MotifResistant_dbProne_TEST &&  mv mast.xml mast_RP_Split_%d_TEST.xml &&  cp mast_RP_Split_%d_TEST.xml	../../../../../MAST && cd /home/data/dicleyalcin/%s_%s/MEME/TrainP/TrainP_Split_%d/meme_out_resistants_split_%d/MAST_outputs_MotifResistant_dbResistant_TEST &&  mv mast.xml mast_RR_Split_%d_TEST.xml && cp mast_RR_Split_%d_TEST.xml ../../../../../MAST"	% (dataset, dataset_property, (i+1),(i+1), (i+1),(i+1),	dataset, dataset_property, (i+1),(i+1),	(i+1),(i+1),dataset, dataset_property, (i+1),(i+1),	(i+1),(i+1),dataset, dataset_property, (i+1),(i+1),	(i+1),(i+1),dataset, dataset_property, (i+1),(i+1),	(i+1),(i+1),dataset, dataset_property, (i+1),(i+1),	(i+1),(i+1),dataset, dataset_property, (i+1),(i+1),	(i+1),(i+1),dataset, dataset_property, (i+1),(i+1), (i+1),(i+1)))
	"""

	#------------------------------------------------------------------------------------------------------------------
	# # # # # # # # # # # # # # #------ unique MAST xml files are created and copied in MAST (overfit case) ----- # # # # # # # # # 
	#------------------------------------------------------------------------------------------------------------------
	

	#return os.system("PATH=$HOME/meme/bin:$PATH && cd /home/dyalcin/data/%s_%s/MEME-overfit/MAST_outputs_MotifProne_dbProne_split_%d/ &&  cp mast.xml mast_PP_Split_%d.xml &&  cp mast_PP_Split_%d.xml /home/dyalcin/data/%s_%s/MAST-overfit && cd /home/dyalcin/data/%s_%s/MEME-overfit/MAST_outputs_MotifProne_dbResistant_split_%d/ && cp mast.xml mast_PR_Split_%d.xml && cp mast_PR_Split_%d.xml /home/dyalcin/data/%s_%s/MAST-overfit && cd /home/dyalcin/data/%s_%s/MEME-overfit/MAST_outputs_MotifResistant_dbProne_split_%d/ && cp mast.xml mast_RP_Split_%d.xml && cp mast_RP_Split_%d.xml /home/dyalcin/data/%s_%s/MAST-overfit	&& cd /home/dyalcin/data/%s_%s/MEME-overfit/MAST_outputs_MotifResistant_dbResistant_split_%d/ && cp mast.xml mast_RR_Split_%d.xml && cp mast_RR_Split_%d.xml /home/dyalcin/data/%s_%s/MAST-overfit	&& 	cd /home/dyalcin/data/%s_%s/MEME-overfit/MAST_outputs_MotifProne_dbProne_TEST_split_%d/ &&  cp mast.xml mast_PP_Split_%d_TEST.xml &&  cp mast_PP_Split_%d_TEST.xml /home/dyalcin/data/%s_%s/MAST-overfit && cd /home/dyalcin/data/%s_%s/MEME-overfit/MAST_outputs_MotifProne_dbResistant_TEST_split_%d/ && cp mast.xml mast_PR_Split_%d_TEST.xml && cp mast_PR_Split_%d_TEST.xml /home/dyalcin/data/%s_%s/MAST-overfit && cd /home/dyalcin/data/%s_%s/MEME-overfit/MAST_outputs_MotifResistant_dbProne_TEST_split_%d/ && cp mast.xml mast_RP_Split_%d_TEST.xml && cp mast_RP_Split_%d_TEST.xml /home/dyalcin/data/%s_%s/MAST-overfit && 	cd /home/dyalcin/data/%s_%s/MEME-overfit/MAST_outputs_MotifResistant_dbResistant_TEST_split_%d/ && cp mast.xml mast_RR_Split_%d_TEST.xml && cp mast_RR_Split_%d_TEST.xml /home/dyalcin/data/%s_%s/MAST-overfit" % (dataset,dataset_property,(i+1),(i+1),(i+1),dataset,dataset_property, dataset,dataset_property,(i+1),(i+1),(i+1),dataset,dataset_property, dataset,dataset_property,(i+1),(i+1),(i+1),dataset,dataset_property, dataset,dataset_property,(i+1),(i+1),(i+1),dataset,dataset_property, dataset,dataset_property,(i+1),(i+1),(i+1),dataset,dataset_property, dataset,dataset_property,(i+1),(i+1),(i+1),dataset,dataset_property, dataset,dataset_property,(i+1),(i+1),(i+1),dataset,dataset_property, dataset,dataset_property,(i+1),(i+1),(i+1),dataset,dataset_property))
	



	#------------------------------------------------------------------------------------------------------------------
	# # # # # # # # # # # # # # # # # ------------ FIMO for TFBS features (UniPROBE) ------------ # # # # # #------------------------------------------------------------------------------------------------------------------
	
	#return os.system("PATH=$HOME/meme/bin:$PATH && cd /home/dyalcin/data/%s_%s/MEME/TrainP/TrainP_Split_%d && fimo --o FIMO_outputs_%s --verbosity 1 /home/dyalcin/data/%s_%s/Motif_databases_MEME/TFBSshape/TFBSshape_%s.meme TrainP_Split_%d.fa && cd /home/dyalcin/data/%s_%s/MEME/TrainR/TrainR_Split_%d &&  fimo --o FIMO_outputs_%s --verbosity 1 /home/dyalcin/data/%s_%s/Motif_databases_MEME/TFBSshape/TFBSshape_%s.meme TrainR_Split_%d.fa && cd /home/dyalcin/data/%s_%s/MEME/TestP/TestP_Split_%d &&  fimo --o FIMO_outputs_%s --verbosity 1 /home/dyalcin/data/%s_%s/Motif_databases_MEME/TFBSshape/TFBSshape_%s.meme TestP_Split_%d.fa && cd /home/dyalcin/data/%s_%s/MEME/TestR/TestR_Split_%d &&  fimo --o FIMO_outputs_%s --verbosity 1 /home/dyalcin/data/%s_%s/Motif_databases_MEME/TFBSshape/TFBSshape_%s.meme TestR_Split_%d.fa" % (dataset,dataset_property,(i+1),TFBS_database, dataset,dataset_property,TFBS_database,(i+1), dataset,dataset_property,(i+1),TFBS_database, dataset,dataset_property,TFBS_database,(i+1), dataset,dataset_property,(i+1),TFBS_database, dataset,dataset_property,TFBS_database,(i+1), dataset,dataset_property,(i+1),TFBS_database, dataset,dataset_property,TFBS_database,(i+1)))
	

	#------------------------------------------------------------------------------------------------------------------
	# # # # # # ------  AFTER FIMO, convert cisml.xml, and move them in TFBS_%s and TFBS_%s_TEST folder. # # # # # # # #     
	#------------------------------------------------------------------------------------------------------------------
	

	
	#return os.system("PATH=$HOME/meme/bin:$PATH && cd /home/dyalcin/data/%s_%s/MEME/TrainP/TrainP_Split_%d/FIMO_outputs_%s && cp cisml.xml cisml_pdb_split_%d.xml && cp cisml_pdb_split_%d.xml /home/dyalcin/data/%s_%s/TFBS_%s/ && cd /home/dyalcin/data/%s_%s/MEME/TrainR/TrainR_Split_%d/FIMO_outputs_%s && cp cisml.xml cisml_rdb_split_%d.xml &&  cp cisml_rdb_split_%d.xml /home/dyalcin/data/%s_%s/TFBS_%s/ && cd /home/dyalcin/data/%s_%s/MEME/TestP/TestP_Split_%d/FIMO_outputs_%s && cp cisml.xml cisml_test_pdb_split_%d.xml && cp cisml_test_pdb_split_%d.xml /home/dyalcin/data/%s_%s/TFBS_%s_TEST/ && cd /home/dyalcin/data/%s_%s/MEME/TestR/TestR_Split_%d/FIMO_outputs_%s && cp cisml.xml cisml_test_rdb_split_%d.xml && cp cisml_test_rdb_split_%d.xml /home/dyalcin/data/%s_%s/TFBS_%s_TEST/" % (dataset,dataset_property,(i+1), TFBS_database,(i+1),(i+1),dataset,dataset_property,TFBS_database,dataset,dataset_property,(i+1),TFBS_database,(i+1),(i+1),dataset, dataset_property,TFBS_database,dataset,dataset_property,(i+1),TFBS_database,(i+1),(i+1),dataset,dataset_property,TFBS_database, dataset,dataset_property,(i+1),TFBS_database,(i+1),(i+1),dataset,dataset_property,TFBS_database))
	

	
	




Parallel(n_jobs=num_cores)(delayed(processInput)(i) for i in inputs)
