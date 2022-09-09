% This template can be used to run q-resolved reconstructions.
% First prepare a STEP2 template for the reconstruction with the staps that
% you want followed and the initial guess that is appropriate.
SASTT_template = 'STEP2_SASTT_optimization_SH_tomo1_reg3_MGS';
indices_q = [100:110];

for indsq = indices_q
    run(SASTT_template)
end



