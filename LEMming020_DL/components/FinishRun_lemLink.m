% FinishRun_lemLink.m - post-run plotting and cleanup
% Save landscape matrix and specified variables to the run's subfolder
% -------------------------------------------------------------------------

disp('Saving final state... ')

% Save the landscape matrix
save(['./' run_filename '/' run_name '_EndState.mat']);

disp('Saved.')

if SAVEMODE == 1 || SAVEMODE == 2

    % Prompt to run MakeLEMmingMov.m
    makemovie = input('Do you want to make a movie of this run? Y/N [Y]: ', 's');
    if isempty(makemovie) || makemovie == 'y';
        makemovie = 'Y';
    end

else
    disp('Done!')
end

if makemovie == 'Y'
    MakeLEMmingMov
else
    disp('Done!')
end