
import fasta
import bowtie_search as bts
import time

hot = fasta.Fasta("HOTAIR.txt") 
tstart = time.time()
hits = bts.align(hot,16,'humanPseudo')

print "%f seconds" % (time.time() - tstart )
