function pdcchCandidates = generate_candidate(candidate_len,total_len)
    pdcchCandidates = [1,candidate_len];
    for i = 73:72:total_len - candidate_len
        tmpcandidate = [i,i+candidate_len-1];
        pdcchCandidates = [pdcchCandidates; tmpcandidate];
    end
end