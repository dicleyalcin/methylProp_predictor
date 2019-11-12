import xml.etree.ElementTree
from xml.dom import minidom
import xlsxwriter
import pdb




for k in range(500):

	xmldoc_P = xml.etree.ElementTree.parse('meme_trainP_split_%d.xml' % (k+1)).getroot()
	xmldoc_R = xml.etree.ElementTree.parse('meme_trainR_split_%d.xml' % (k+1)).getroot()






	workbook = xlsxwriter.Workbook('SigMotif_Split_%d.xlsx' %(k+1))

	worksheet1 = workbook.add_worksheet('ProneMotifs')
	worksheet2 = workbook.add_worksheet('ResistantMotifs')



	# MOTIF NAMES (columns)

	for i in range(30):
		worksheet1.write('A%d' % (i+2), 'PM%d' %(i+1))
		worksheet2.write('A%d' % (i+2), 'RM%d' %(i+1))




	evalue_P=[]
	for motif in xmldoc_P.iter('motif'):
		a = motif.attrib
		evalue = a.get('e_value')
		evalue_P.append(evalue)
	for j in range(30):
		worksheet1.write('B%d' % (j+2), evalue_P[j])


	evalue_R=[]
	for motif in xmldoc_R.iter('motif'):
		a = motif.attrib
		evalue = a.get('e_value')
		evalue_R.append(evalue)
	for j in range(30):
		worksheet2.write('B%d' % (j+2), evalue_R[j])

	workbook.close()
