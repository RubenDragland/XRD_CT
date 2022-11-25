function [finish] = run_dummies()
    dummies = ["Difficult_periodic_filter3_17cube_8off_0align" "Dummy_periodic_filter1_17cube_8off_0align" "Dummy_periodic_filter1_27cube_8off_0align" "Dummy_periodic_filter3_17cube_18off_0align"];
    filters = [3 1 1 1];
    dir = cd;

    for i = 1:length(dummies)
        Run_SH(dir, dummies(i), 0, filters(i));
        Run_SH(dir, dummies(i), 1, filters(i));   
    end
end 