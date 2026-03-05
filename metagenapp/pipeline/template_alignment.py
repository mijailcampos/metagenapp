def align_to_template(query_fasta, template_aln_fasta, output_fasta,
                      k=9, top_candidates=30,
                      min_id=0.75, min_cov=0.75,
                      threads=8):

    template = load_template(template_aln_fasta)
    # template.raw_seqs[], template.aln_seqs[], template.raw_to_aln_maps[]
    kmer_index = build_kmer_index(template.raw_seqs, k=k, max_hits_per_kmer=500)

    with open(output_fasta, "w") as out:
        for q_id, q_seq in stream_fasta(query_fasta):
            cands = pick_candidates(q_seq, kmer_index, k=k, top=top_candidates)

            best = None
            for tid in cands:
                aln = align(q_seq, template.raw_seqs[tid])  # edlib ideal
                score = compute_score(aln)
                best = max(best, (score, tid, aln))

            if best is None or best.score < thresholds:
                continue  # o escribir a un archivo unmapped

            aligned_query = project_to_template(
                q_seq,
                best.aln,
                template.raw_to_aln_maps[best.tid],
                template.aln_length
            )

            out.write(f">{q_id}\n{aligned_query}\n")