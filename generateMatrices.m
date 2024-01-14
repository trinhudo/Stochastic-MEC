function matrices = generateMatrices()
rows = 3;
cols = 2;

% Initialize matrices cell array
matrices = cell(1, nchoosek(cols, 1)^rows);

% Generate all permutations of length `cols` for each row
permutations = perms(1:cols);

% Generate all combinations of length `rows` from the permutations
combinations = combvec(permutations(:, 1), permutations(:, 2));

% Reshape the combinations into 3x2 matrices
matrices = reshape(combinations, rows, cols, []);

% Convert to a cell array for better display
matrices = num2cell(matrices, [1, 2]);

% Display the matrices
disp('All possible matrices:');
for i = 1:length(matrices)
    disp(matrices{i});
end
end
