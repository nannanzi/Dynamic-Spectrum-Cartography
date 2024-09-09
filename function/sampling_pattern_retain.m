function SampleIndex = sampling_pattern_retain(I,J,p,ratio)
%input: I (length); J(width); p(sampling_ratio)
    AvaiIndex = 1:I*J;
    SampPerIt = round(length(AvaiIndex)*p);
    
    SampleNew = randperm(numel(AvaiIndex),SampPerIt);
    SampleIndex{1} = SampleNew;
    AvaiIndex(SampleNew) = [];
    tt = 2;
    while numel(AvaiIndex) > 0
        
        IndexRetain = randperm(SampPerIt,round(SampPerIt*ratio));
        SampleRetain = SampleIndex{tt-1}(IndexRetain);
        
        if round(SampPerIt*(1 - ratio)) > numel(AvaiIndex)
            SampleRetain = SampleIndex{tt-1};
            IndexNew = 1:numel(AvaiIndex);
            SampleNew = AvaiIndex;
        else
            IndexNew = randperm(numel(AvaiIndex),round(SampPerIt*(1 - ratio)));
            SampleNew = AvaiIndex(IndexNew);
        end
        
        AvaiIndex(IndexNew) = [];
        
        SampleIndex{tt} = [SampleRetain,SampleNew];
        tt = tt + 1;
    end
    

end