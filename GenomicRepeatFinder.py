import re
from Bio import SeqIO
import sys
import string

BIN_SIZE = 1000000
FORWARD_SEQ = 'TTTAGGG'
MIN_REPEATS = 3
FASTA_FILE = []

## To use sliding window opposed to fixed window, alter Window & Step sizes
WINDOW_SIZE = BIN_SIZE
STEP_SIZE = WINDOW_SIZE

def prompt_user():
	global FASTA_FILE, WINDOW_SIZE, STEP_SIZE, MIN_REPEATS, FORWARD_SEQ
	if len(sys.argv) != 2:
		FASTA_FILE = input("Enter the path to the FASTA file: ")
	else:
		FASTA_FILE = sys.argv[1] or input("Enter the path to the FASTA file: ")
	FORWARD_SEQ = input(f"Enter the forward pattern (default is {FORWARD_SEQ}):").strip().upper() or FORWARD_SEQ	
	try:
		MIN_REPEATS = int(input(f"Enter minimum repeats (default is {MIN_REPEATS}): ").strip() or MIN_REPEATS)
	except ValueError:
		print("Invalid input. Using default minimum repeats.")
##	WINDOW_SIZE = (len(FORWARD_SEQ)) * MIN_REPEATS
	try:
		WINDOW_SIZE = int(input(f"Enter window size (default is {WINDOW_SIZE}): ").strip() or WINDOW_SIZE)
	except ValueError:
		print("Invalid input. Using default window size.")
##	STEP_SIZE = WINDOW_SIZE - (len(FORWARD_SEQ) * MIN_REPEATS - 1)
	print(f"Step Size set to {STEP_SIZE}")
def process_all_records(fasta_records, FORWARD_SEQ, MIN_REPEATS, STEP_SIZE, WINDOW_SIZE, define_pattern, circle_match):
    for fasta_record in fasta_records:
        fasta = str(fasta_record.seq).upper()
        fasta_length = len(fasta)
        reverse_seq = complement(FORWARD_SEQ)
        forward_pattern = '|'.join(define_pattern(FORWARD_SEQ))
        reverse_pattern = '|'.join(define_pattern(reverse_seq))

        NUM_BINS = (fasta_length + BIN_SIZE - 1) // BIN_SIZE

        # Initialize observations with 0 for forward_hits and reverse_hits
        observations = {bin_num: {'forward_hits': 0, 'reverse_hits': 0} for bin_num in range(NUM_BINS)}

        with open(f"{fasta_record.id}_R{MIN_REPEATS}.txt", 'w') as matches_out:
            matches_out.write('id\tgroup\tobservation\tvalue\n')

            for i in range(0, fasta_length + 1, STEP_SIZE):
                window = fasta[i:i + WINDOW_SIZE]
                fhits, fmatches = circle_match(forward_pattern, window, i, MIN_REPEATS)
                rhits, rmatches = circle_match(reverse_pattern, window, i, MIN_REPEATS)

                for match_start, match_seq in fmatches.items():
                    bin_number = get_bin(match_start)
                    observations[bin_number]['forward_hits'] += int(len(match_seq) / len(FORWARD_SEQ))

                for match_start, match_seq in rmatches.items():
                    bin_number = get_bin(match_start)
                    observations[bin_number]['reverse_hits'] -= int(len(match_seq) / len(reverse_seq))

            for bin_number in range(NUM_BINS):
                start_position = bin_number * BIN_SIZE
                end_position = min(start_position + BIN_SIZE, fasta_length)
                matches_out.write(f'{start_position}\t{fasta_record.id}\tforward_hits\t{observations[bin_number]["forward_hits"]}\n')
                matches_out.write(f'{start_position}\t{fasta_record.id}\treverse_hits\t{observations[bin_number]["reverse_hits"]}\n')

        print(f"finished processing {fasta_record.id}")
def complement(sequence):
	code = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
	complement = []
	for base in reversed(sequence):
		if base in code:
			complement.append(code[base])
		else:
			print(f"Warning: Unable to generate complement base to '{base}' in the sequence.")
			return None
	return ''.join(complement)
def define_pattern(pattern):
	final_pattern = []
	for ix in range(len(pattern)):
		new_pattern = pattern[ix:] + "(" + pattern + ")" + "{" + str(MIN_REPEATS - 1) + ",}" + pattern[0:ix]
		final_pattern.append(new_pattern)
	return final_pattern
def circle_match(pattern, fasta, pos, MIN_REPEATS):
	pattern_length = len(pattern.split('|')[0].split('(')[0])
	hits = 0
	matches = {}
	for match in re.finditer(pattern, fasta):
		match_length = len(match.group(0))
		if match_length >= pattern_length * MIN_REPEATS:
			repeat_count = match_length / pattern_length
			hits += repeat_count
			matches[match.start() + pos] = match.group(0)
	return hits, matches
def get_bin(position):
    return position // BIN_SIZE

if __name__ == "__main__":
    prompt_user()
    try:
        fasta_records = list(SeqIO.parse(FASTA_FILE, format="fasta"))
        process_all_records(fasta_records, FORWARD_SEQ, MIN_REPEATS, STEP_SIZE, WINDOW_SIZE, define_pattern, circle_match)
    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)
