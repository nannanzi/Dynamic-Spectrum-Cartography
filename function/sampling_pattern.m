% Function: sampling pattern. Caveat: guarantee that each column/row has at least one sample
function SampleIndex = sampling_pattern(I,J,p)
%input: I (length); J(width); p(sampling_ratio)
    basic_sampling_ratio = 2/J;
    SampleIndex = [];
    sampled_row_per_column = ceil(I*(p-basic_sampling_ratio));
    
    for ii = 1:I
        row_idx = randperm(J,2);
        sample_idx_ii= ii + (row_idx-1)*I;
        SampleIndex = [SampleIndex,sample_idx_ii];
    end
    
    
    for jj = 1:J
        col_idx = randperm(I,sampled_row_per_column);
        sample_idx_jj = col_idx + (jj-1)*I;
        SampleIndex = [SampleIndex,sample_idx_jj];
    end
    SampleIndex = SampleIndex';
end