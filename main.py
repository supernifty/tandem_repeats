
import argparse
import collections
import logging
import sys

def find_run(data, kmer, start, repeat):
  result = 0
  while kmer == find_kmer_at_pos(data, start, repeat):
    result += repeat
    start += repeat
    if start >= len(data) - repeat:
      return result
  return result

def find(fh, repeat, minlen):
  '''
    find tandem repeats of at least minlen
  '''
  #sys.stdout.write('{},{},{},{}\n'.format("chrom", "pos", "length", "kmer"))
  last_log = 0
  chrom_data = []
  for line in sys.stdin:
    if line.startswith('>'):
      if len(chrom_data) > 0:
        find_tandems(chrom, ''.join(chrom_data), repeat, minlen)
      # restart
      chrom = line.strip('\n')[1:]
      logging.info('reading %s...', chrom)
      last = ''
      chrom_data = []
    else:
      line = line.strip('\n').upper()
      combined = last + line
      chrom_data.append(line)

  if len(chrom_data) > 0:
    find_tandems(chrom, ''.join(chrom_data), repeat, minlen)

'''
  what kmer is at pos?
  check each kmer and see if it is in the set
  this isn't efficient.
'''
def find_kmer_at_pos(data, pos, repeat):
  kmer = data[pos:pos+repeat]
  #logging.debug('%s %i %i', kmer, pos, len(data))
  if 'N' in kmer:
    return None
  if repeat > 1 and all([kmer[x] == kmer[0] for x in range(1, repeat)]): # mono
    return None
  if repeat > 2 and repeat % 2 == 0 and kmer[:repeat / 2] == kmer[repeat / 2:]: # abab
    return None
  return kmer

def find_tandems(chrom, chrom_data, repeat, minlen):
  # now find tandem repeats
  logging.info('finding tandem repeats for %s at %i positions', chrom, len(chrom_data))
  max_len = 0
  max_len_kmer = collections.defaultdict(int)
  run_to = 0
  last_log = 0
  for pos in range(len(chrom_data) - repeat + 1): #sorted(kmers[kmer]):
    kmer = find_kmer_at_pos(chrom_data, pos, repeat)
    if kmer is None:
      continue
    if pos < run_to: 
      continue
    run_len = find_run(chrom_data, kmer, pos, repeat)
    #logging.debug('run_len at %i is %i with %s', pos, run_len, kmer)
    max_len = max(max_len, run_len)
    max_len_kmer[kmer] = max(max_len_kmer[kmer], run_len)
    if run_len >= minlen:
      #sys.stdout.write('{},{},{},{}\n'.format(chrom, pos, run_len, kmer))
      sys.stdout.write('{}\t{}\t{}\trepeat={};length={}\n'.format(chrom, pos, pos + run_len, kmer, run_len))
      run_to = pos + run_len
    if pos - last_log > 10e6:
      logging.debug('checked %s:%i...', chrom, pos)
      last_log = pos
  logging.debug('max run found for kmers: %s', ', '.join(['{}: {}'.format(kmer, max_len_kmer[kmer]) for kmer in sorted(max_len_kmer.keys())]))
  logging.debug('max run found for %s: %i', chrom, max_len)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Find tandem repeats')
  parser.add_argument('--repeat', type=int, default=1, help='find tandem repeats of this repeat size')
  parser.add_argument('--min', type=int, default=4, help='only report tandem repeats of this length')
  args = parser.parse_args()
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG) 
  find(sys.stdin, args.repeat, args.min)
