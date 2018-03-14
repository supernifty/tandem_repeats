
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
        find_tandems(chrom, kmers, repeat, minlen)
      # restart
      chrom = line.strip('\n')[1:]
      logging.info('building dictionary for %s', chrom)
      pos = 0
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
      pos += len(line)
      last = line[-repeat:]
      if pos - last_log > 10e6:
        logging.debug('read %s:%i...', chrom, pos)
        last_log = pos

  if pos > 0:
    find_tandems(chrom, kmers, repeat, minlen)

def find_tandems(chrom, kmers, repeat, minlen):
  # now find tandem repeats
  logging.info('finding tandem repeats for %s', chrom)
  max_len = 0
  for kmer in kmers:
    logging.debug('finding tandem repeats for %s. evaluating %s with %i locations', chrom, kmer, len(kmers[kmer]))
    max_len_kmer = 0
    run_to = 0
    for pos in sorted(kmers[kmer]):
      if pos < run_to: 
        continue
      run_len = find_run(kmers[kmer], pos, repeat)
      max_len = max(max_len, run_len)
      max_len_kmer = max(max_len_kmer, run_len)
      if run_len >= minlen:
        sys.stdout.write('{},{},{},{}\n'.format(chrom, pos, run_len, kmer))
      run_to = pos + run_len
    logging.debug('max run found for %s.%s: %i', chrom, kmer, max_len_kmer)
  logging.debug('max run found for %s: %i', chrom, max_len)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Compare BAMs')
  parser.add_argument('--repeat', type=int, default=1, help='find tandem repeats of this repeat size')
  parser.add_argument('--min', type=int, default=10, help='only report tandem repeats of this length')
  args = parser.parse_args()
  logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s', level=logging.DEBUG) 
  find(sys.stdin, args.repeat, args.min)
