import io
import os
from repDNA.util import get_data
from repDNA.nac import Kmer
import numpy as np
import itertools


bases=['A','T','G','C']

for i in range(9):
	header = [''.join(p) for p in itertools.product(bases, repeat=i+2)]
	kmer=Kmer(i+2,normalize=False,upto=False)
	pos_vec = kmer.make_kmer_vec(open('/home/dyalcin/data/raw_datas/rawdata_estecio/prones.fa'))
	neg_vec = kmer.make_kmer_vec(open('/home/dyalcin/data/raw_datas/rawdata_estecio/resistants.fa'))

	with io.open("PronesK=" + str(i+2) + ".txt", 'w+') as f:
		f.write(unicode(header))
		f.write(unicode(pos_vec))
		f.close()

	with io.open("ResistantsK=" + str(i+2) + ".txt", 'w+',) as f2:
		f2.write(unicode(neg_vec))
		f2.close()
