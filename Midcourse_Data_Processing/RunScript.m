% % % INSTRUCTIONS: The first argument off the batch() function should be the matlab script to run.
% % % If multiple CPUs are used, enter the count as the value of 'Pool'.

c=parcluster;
c.NumWorkers = 6;
saveProfile(c)

% % Use for multiple CPUs
my_job = batch ( 'ClusterAlgorithmParallelTermination', 'Profile', 'local' , 'Pool', 5);

% % % % Use for single CPU
%my_job = batch ( 'ClusterAlgorithmParallelTermination', 'Profile', 'local');

wait ( my_job );

diary ( my_job );

% % % ==================================================================
% % % Extract error output from each worker (CPU)
% % % ==================================================================

fid = fopen('error_log.txt', 'w');

for i = 1:length(my_job.Tasks)
    if ~isempty(my_job.Tasks(i).Error)
        error_details = getReport(my_job.Tasks(i).Error);
        % write to file
        fprintf(fid, 'Error in task %d:\n%s\n\n', i, error_details);
    end
end

fclose(fid);


delete ( my_job );