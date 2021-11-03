import os
import codecs
import collections
import argparse
from collections import defaultdict
from itertools import chain
from collections import Counter
from operator import itemgetter
import json
import sys

GAP_ERRORS = True
SOURCE_ERRORS = True


def read_file(file_path, alignments=False):
    with codecs.open(file_path, 'r', "utf-8") as fid:
        lines = [line.rstrip().split() for line in fid.readlines()]
    if alignments:
        alignments = []
        for sent in lines:
            alignments.append([
                tuple([
                    int(side) if side else None for side in word.split('-')
                ])
                for word in sent
            ])
        return alignments
    else:
        return lines


def check_out_of_bounds(tokens, alignments, source=True):
    """
    Checks if alignment indices our out of bounds to respect to tokens. This
    can happend as we generated the alignments with tercom but use the original
    raw files (encoding problems may lead to tokens disspearing)
    """
    assert len(tokens) == len(alignments), "Number of sentences does not macth"
    for sent_index in range(len(tokens)):

        length = len(tokens[sent_index])

        # Adaption for Python3
        # if source:
        #     max_index = max([x[0] for x in alignments[sent_index]])
        # else:
        #     max_index = max([x[1] for x in alignments[sent_index]])

        if source:
            max_index = max([x[0] if x[0] is not None else -1 for x in alignments[sent_index]])
        else:
            max_index = max([x[1] if x[1] is not None else -1 for x in alignments[sent_index]])
        if max_index >= length:
            print("Sentence Index: %d" % sent_index)
            print(tokens[sent_index])
            print(alignments[sent_index])
            raise Exception(
                "Tercom alignments and original tokens do not match."
                "Likely an enconding problem"
            )


def parse_arguments(sys_argv):
    parser = argparse.ArgumentParser(
        prog='New word-level tags version 0.0.1'
    )
    # Arguments defining a time slice
    parser.add_argument(
        '--in-source-tokens',
        help='One sentence per line',
        required=True,
        type=str
    )
    parser.add_argument(
        '--in-mt-tokens',
        help='One sentence per line',
        required=True,
        type=str
    )
    parser.add_argument(
        '--in-pe-tokens',
        help='One sentence per line',
        required=True,
        type=str
    )
    parser.add_argument(
        '--in-source-pe-alignments',
        help='Fast align format',
        required=True,
        type=str
    )
    parser.add_argument(
        '--in-pe-mt-alignments',
        help='Fast align format, deletions insertions as empty index',
        required=True,
        type=str
    )

    # OUTPUTS
    parser.add_argument(
        '--out-source-tags',
        help='Source OK/BAD tags per sentence',
        required=True,
        type=str
    )
    parser.add_argument(
        '--out-target-tags',
        help='Target OK/BAD tags per sentence',
        required=True,
        type=str
    )

    # new argument indicating output file of SRC-MT alignments
    parser.add_argument(
        '--out-source-mt-word-alignments',
        help='Source-MT Word alignment per sentence.',
        default=None,
        type=str
    )
    parser.add_argument(
        '--out-source-mt-gap-alignments',
        help='Source-MT Gap alignment per sentence.',
        default=None,
        type=str
    )

    parser.add_argument(
        '--fluency-rule',
        help='Rules used to determine source tags',
        choices=['ignore-shift-set', 'normal', 'missing-only'],
        type=str
    )
    args = parser.parse_args(sys_argv)

    return args


def read_data(args):
    source_tokens = read_file(args.in_source_tokens)
    mt_tokens = read_file(args.in_mt_tokens)
    pe_tokens = read_file(args.in_pe_tokens)
    src_pe_alignments = read_file(args.in_source_pe_alignments, alignments=True)
    pe_mt_alignments = read_file(args.in_pe_mt_alignments, alignments=True)

    # Sanity Checks
    # Number of sentences matches
    num_sentences = len(source_tokens)
    assert len(mt_tokens) == num_sentences, \
        "Number of sentences does not match"
    assert len(pe_tokens) == num_sentences, \
        "Number of sentences does not match"
    assert len(src_pe_alignments) == num_sentences, \
        "Number of sentences does not match"
    assert len(pe_mt_alignments) == num_sentences, \
        "Number of sentences does not match"
    # fast_align alignments out-of-bounds
    check_out_of_bounds(source_tokens, src_pe_alignments, source=True)
    check_out_of_bounds(pe_tokens, src_pe_alignments, source=False)
    # tercom alignments out-of-bounds
    check_out_of_bounds(mt_tokens, pe_mt_alignments, source=False)
    check_out_of_bounds(pe_tokens, pe_mt_alignments, source=True)

    # Reorganize source-target alignments as a dict
    pe2source = []
    for sent in src_pe_alignments:
        pe2source_sent = defaultdict(list)
        for src_idx, pe_idx in sent:
            pe2source_sent[pe_idx].append(src_idx)
        pe2source.append(pe2source_sent)

    return (
        source_tokens, mt_tokens, pe_tokens, pe2source, pe_mt_alignments
    )


def get_quality_tags(mt_tokens, pe_tokens, pe_mt_alignments, pe2source,
                     fluency_rule=None):
    # Word + Gap Tags
    target_tags = []
    source_tags = []
    error_detail = []

    # SRC-MT Alignment
    src_mt_word_alignments = []
    src_mt_gap_alignments = []

    for sentence_index in range(len(mt_tokens)):

        # Variables for this sentence
        sent_tags = []
        sent_deletion_indices = []
        source_sentence_bad_indices = set()
        error_detail_sent = []
        mt_position = 0
        source_mt_word_alignment = collections.defaultdict(set)
        source_mt_gap_alignment = collections.defaultdict(set)

        # Loop over alignments. This has the length of the edit-distance aligned
        # sequences.
        prev_mt_idx = None
        for pe_idx, mt_idx in pe_mt_alignments[sentence_index]:

            if mt_idx is None:

                # Deleted word error (need to store for later)
                sent_deletion_indices.append(mt_position - 1)

                if fluency_rule == 'normal' or fluency_rule == "missing-only":

                    source_positions = pe2source[sentence_index][pe_idx]
                    source_sentence_bad_indices |= set(source_positions)
                    error_type = 'deletion'

                    # yield Source-GAP Alignment
                    for source_pos in source_positions:
                        source_mt_gap_alignment[source_pos].add(mt_position)

                elif fluency_rule == 'ignore-shift-set':

                    # RULE: If word exists elsewhere in the sentence do not
                    # propagate error to the source.
                    if (
                            pe_tokens[sentence_index][pe_idx] not in
                            mt_tokens[sentence_index]
                    ):
                        source_positions = pe2source[sentence_index][pe_idx]
                        source_sentence_bad_indices |= set(source_positions)
                        error_type = 'deletion'

                        for source_pos in source_positions:
                            source_mt_gap_alignment[source_pos].add(mt_position)
                    else:
                        source_positions = None
                        error_type = 'deletion (shift)'

                else:
                    raise Exception("Uknown rule %s" % fluency_rule)

                # Store error detail
                error_detail_sent.append({
                    'type': error_type,
                    'gap_position': mt_position - 1,
                    'target_position': mt_idx,
                    'source_positions': source_positions,
                })

            elif pe_idx is None:

                # Insertion error
                if mt_idx == prev_mt_idx:
                    if sent_tags[-1] == 'OK':
                        sent_tags[-1] = 'BAD'
                else:
                    sent_tags.append('BAD')
                    mt_position += 1

                # Store error detail
                error_detail_sent.append({
                    'type': 'insertion',
                    'target_position': mt_idx,
                    'source_positions': None,
                })

            elif (
                    mt_tokens[sentence_index][mt_idx] !=
                    pe_tokens[sentence_index][pe_idx]
            ):

                # Substitution error
                if mt_idx == prev_mt_idx:
                    if sent_tags[-1] == 'OK':
                        sent_tags[-1] = 'BAD'
                else:
                    sent_tags.append('BAD')
                    mt_position += 1

                # Aligned words in the source are BAD
                # RULE: If word exists elsewhere in the sentence do not
                # propagate error to the source.
                if fluency_rule == 'normal':

                    source_positions = pe2source[sentence_index][pe_idx]
                    source_sentence_bad_indices |= set(source_positions)
                    error_type = 'substitution'

                    for source_pos in source_positions:
                        source_mt_word_alignment[source_pos].add(mt_idx)

                elif fluency_rule == 'ignore-shift-set':

                    # RULE: If word exists elsewhere in the sentence do not
                    # propagate error to the source.
                    if (
                            pe_tokens[sentence_index][pe_idx] not in
                            mt_tokens[sentence_index]
                    ):
                        source_positions = pe2source[sentence_index][pe_idx]
                        source_sentence_bad_indices |= set(source_positions)
                        error_type = 'substitution'

                        for source_pos in source_positions:
                            source_mt_word_alignment[source_pos].add(mt_idx)
                    else:
                        source_positions = None
                        error_type = 'substitution (shift)'

                elif fluency_rule == 'missing-only':

                    source_positions = None
                    error_type = 'substitution'

                else:
                    raise Exception("Uknown rule %s" % fluency_rule)

                # Store error detail
                error_detail_sent.append({
                    'type': error_type,
                    'target_position': mt_idx,
                    'source_positions': source_positions,
                })

            else:

                # OK
                if mt_idx == prev_mt_idx:
                    continue
                sent_tags.append('OK')
                mt_position += 1

                source_positions = pe2source[sentence_index][pe_idx]
                for source_pos in source_positions:
                    source_mt_word_alignment[source_pos].add(mt_idx)

            if mt_idx is not None and mt_idx != prev_mt_idx:
                prev_mt_idx = mt_idx

        # Insert deletion errors as gaps
        if GAP_ERRORS:
            word_and_gaps_tags = []

            # Add starting OK/BAD
            if -1 in sent_deletion_indices:
                word_and_gaps_tags.append('BAD')
            else:
                word_and_gaps_tags.append('OK')

            # Add rest of OK/BADs
            for index, tag in enumerate(sent_tags):
                if index in sent_deletion_indices:
                    word_and_gaps_tags.extend([tag, 'BAD'])
                else:
                    word_and_gaps_tags.extend([tag, 'OK'])
            target_tags.append(word_and_gaps_tags)
        else:
            target_tags.append(sent_tags)

        # Convert BAD source indices into indices
        source_sentence_bad_tags = \
            ['OK'] * len(source_tokens[sentence_index])
        for index in list(source_sentence_bad_indices):
            source_sentence_bad_tags[index] = 'BAD'
        source_tags.append(source_sentence_bad_tags)

        #
        error_detail.append(error_detail_sent)

        # yield Source-GAP Alignment of one triplets of data
        src_mt_gap_alignments.append(source_mt_gap_alignment)
        src_mt_word_alignments.append(source_mt_word_alignment)

    # Basic sanity checks
    assert all(
        [len(aa) * 2 + 1 == len(bb) for aa, bb in zip(mt_tokens, target_tags)]
    ), "tag creation failed"
    assert all(
        [len(aa) == len(bb) for aa, bb in zip(source_tokens, source_tags)]
    ), "tag creation failed"

    # check the sanity of SRC-MT Gap alignments with generated tags
    for sent_i, src_mt_gap_alignment in enumerate(src_mt_gap_alignments):
        for src_i, tgt_is in src_mt_gap_alignment.items():
            assert source_tags[sent_i][src_i] == 'BAD'
            if source_tags[sent_i][src_i] != 'BAD':
                raise Exception('INS Source word do not have BAD tag.')
            for tgt_i in tgt_is:
                if target_tags[sent_i][tgt_i * 2] != 'BAD':
                    raise Exception('INS MT gap do not have BAD tag.')

    # check the sanity of SRC-MT Word alignment with generated tags
    for sent_i, src_mt_word_alignment in enumerate(src_mt_word_alignments):
        for src_i, tgt_is in src_mt_word_alignment.items():
            for tgt_i in tgt_is:
                src_tag = source_tags[sent_i][src_i]
                tgt_tag = target_tags[sent_i][tgt_i * 2 + 1]
                if src_tag == 'OK' and tgt_tag != 'OK':
                    raise ValueError('Unmatched Tags')
                if tgt_tag == 'BAD' and src_tag != 'BAD':
                    raise ValueError('Unmatched Tags')

    # return source_tags, target_tags, error_detail
    return source_tags, target_tags, error_detail, src_mt_word_alignments, src_mt_gap_alignments


def write_tags(output_file, tags):
    with open(output_file, 'w') as fid:
        for sent_tags in tags:
            tags_line = " ".join(sent_tags)
            fid.write("%s\n" % tags_line)


def write_error_detail(output_file, error_detail):
    with open(output_file, 'w') as fid:
        for error_sent in error_detail:
            fid.write("%s\n" % json.dumps(error_sent))


if __name__ == '__main__':

    # ARGUMENT HANDLING
    args = parse_arguments(sys.argv[1:])

    # READ DATA AND CHECK INTEGRITY
    (
        source_tokens,
        mt_tokens,
        pe_tokens,
        pe2source,
        pe_mt_alignments
    ) = read_data(args)

    # GET TAGS FOR SOURCE AND TARGET
    source_tags, target_tags, error_detail, src_mt_word_alignments, src_mt_gap_alignments = \
        get_quality_tags(
            mt_tokens,
            pe_tokens,
            pe_mt_alignments,
            pe2source,
            fluency_rule=args.fluency_rule
        )

    # Store a more details summary of errors
    error_detail_flat = list(chain.from_iterable(error_detail))
    print(Counter(map(itemgetter('type'), error_detail_flat)))
    dirname = os.path.dirname(args.out_source_tags)
    basename = os.path.basename(args.out_source_tags).split('.')[0]
    error_detail_json = "%s/%s.json" % (dirname, basename)
    write_error_detail(error_detail_json, error_detail)
    print("Wrote %s" % error_detail_json)

    # WRITE DATA
    write_tags(args.out_source_tags, source_tags)
    print("Wrote %s" % args.out_source_tags)
    write_tags(args.out_target_tags, target_tags)
    print("Wrote %s" % args.out_target_tags)

    # Write SRC-MT Word alignments
    if args.out_source_mt_word_alignments:
        with open(args.out_source_mt_word_alignments, 'w') as wf:
            for row_info in src_mt_word_alignments:
                aligns = []
                for k in sorted(row_info):
                    for t_i in row_info[k]: aligns.append(f'{k}-{t_i}')
                wf.write(' '.join(aligns) + '\n')

    # Write SRC-MT Gap alignments
    if args.out_source_mt_gap_alignments:
        with open(args.out_source_mt_gap_alignments, 'w') as wf:
            for row_info in src_mt_gap_alignments:
                aligns = []
                for k in sorted(row_info):
                    for t_i in row_info[k]: aligns.append(f'{k}-{t_i * 2}')
                wf.write(' '.join(aligns) + '\n')