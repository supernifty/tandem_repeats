
import argparse
import collections
import logging
import sys

def find_run(matches, start, repeat):
  result = 0
  while start in matches:
    result += repeat
    start += repeat
  return result

def find(fh, repeat, minlen):
  sys.stdout.write('{},{},{},{}\n'.format("chrom", "pos", "length", "kmer"))
  pos = 0
  last_log = 0
  for line in sys.stdin:
    if line.startswith('>'):
      if pos > 0:
        find_tandems(chrom, kmers, repeat, minlen, max_pos)
      # restart
      chrom = line.strip('\n')[1:]
      logging.info('building dictionary for %s', chrom)
      pos = 0
      max_pos = 0
      last = ''
      kmers = collections.defaultdict(set)
    else:
      line = line.strip('\n').upper()
      combined = last + line
      for offset in range(0, len(combined) - repeat + 1):
        kmer = combined[offset:offset+repeat]
        if 'N' in kmer:
          continue
        if repeat > 1 and all([kmer[x] == kmer[0] for x in range(1, repeat)]):
          continue
        kmers[kmer].add(pos - len(last) + offset)
        max_pos = pos - len(last) + offset
      pos += len(line)
      last = line[-repeat:]
      if pos - last_log > 10e6:
        logging.debug('read %s:%i...', chrom, pos)
        last_log = pos

  if pos > 0:
    find_tandems(chrom, kmers, repeat, minlen, max_pos)

def find_kmer_at_pos(kmers, pos):
  for kmer in kmers:
    if pos in kmers[kmer]:
      return kmer
  return None

def find_tandems(chrom, kmers, repeat, minlen, max_pos):
  # now find tandem repeats
  logging.info('finding tandem repeats for %s at %i positions', chrom, max_pos + 1)
  max_len = 0
  max_len_kmer = collections.defaultdict(int)
  run_to = 0
  last_log = 0
  for pos in range(max_pos + 1): #sorted(kmers[kmer]):
    kmer = find_kmer_at_pos(kmers, pos)
    if kmer is None:
      continue
    if pos < run_to: 
      continue
    run_len = find_run(kmers[kmer], pos, repeat)
    max_len = max(max_len, run_len)
    max_len_kmer[kmer] = max(max_len_kmer[kmer], run_len)
    if run_len >= minlen:
      sys.stdout.write('{},{},{},{}\n'.format(chrom, pos, run_len, kmer))
    run_to = pos + run_len
    if pos - last_log > 10e6:
      logging.debug('checked %s:%i...', chrom, pos)
      last_log = pos
  logging.debug('max run found for kmers: %s', ', '.join(['{}: {}'.format(kmer, max_len_kmer[kmer]) for kmer in kmers]))
  logging.debug('max run found for %s: %i', chrom, max_len)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Compare BAMs')
  parser.add_argument('--repeat', type=int, default=1, help='find tandem repeats of this repeat size')
  parser.add_argument('--min', type=int, default=4, help='only report tandem repeats of this length')
  args = parser.parse_args()
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG) 
  find(sys.stdin, args.repeat, args.min)
