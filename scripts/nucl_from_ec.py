import sys
from Bio import Entrez, SeqIO

Entrez.email = 'jflater@iastate.edu'

# First, find entries that contain the E.C. number
ec_num = sys.argv[1].strip()
# Editing to only search for #, exclude "EC and E.C."
esearch_handle = Entrez.esearch(db='nucleotide', term=ec_num)
entries = Entrez.read(esearch_handle)
esearch_handle.close()

# Second, fetch these entries
efetch_handle = Entrez.efetch(db='nucleotide', id=entries['IdList'], rettype='gb', retmode='xml') 
records = Entrez.parse(efetch_handle)

# Now, we go through the records and look for a feature with name 'EC_number'
for record in records:
      for feature in record['GBSeq_feature-table']:
          for subfeature in feature['GBFeature_quals']:
              if (subfeature['GBQualifier_name'] == 'EC_number'   and
                subfeature['GBQualifier_value'] == ec_num):

                    # If we found it, we extract the seq's start and end
                    accession = record['GBSeq_primary-accession']
                    interval = feature['GBFeature_intervals'][0]
                    interval_start = interval['GBInterval_from']
                    interval_end = interval['GBInterval_to']
                    location = feature['GBFeature_location']
                    if location.startswith('complement'):
                        strand = 2
                    else:
                        strand = 1

                    # Now we fetch the nucleotide sequence
                    handle = Entrez.efetch(db="nucleotide", id=accession,
                                           rettype="fasta", strand=strand,
                                           seq_start = interval_start,
                                           seq_stop = interval_end)
                    seq = SeqIO.read(handle, "fasta")

                    print('>GenBank Accession:{}'.format(accession))
                    print(seq.seq)
efetch_handle.close()
