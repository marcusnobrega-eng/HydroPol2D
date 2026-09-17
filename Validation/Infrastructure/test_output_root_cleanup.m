function tests = test_output_root_cleanup
%TEST_OUTPUT_ROOT_CLEANUP Regression tests for output-directory handling.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
modelRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
functionsDir = fullfile(modelRoot, 'HydroPol2D_Functions');
addpath(functionsDir);
testCase.TestData.FunctionsDir = functionsDir;
end

function teardownOnce(testCase)
rmpath(testCase.TestData.FunctionsDir);
end

function testDeletesReadOnlyNestedResults(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() remove_if_present(root)); %#ok<NASGU>

nested = fullfile(root, 'Modeling_Results', 'Tables_CSV');
mkdir(nested);
oldResult = fullfile(nested, 'old_result.csv');
write_text(oldResult, 'old result');
fileattrib(oldResult, '-w', 'u');

actualRoot = hydropol2d_prepare_output_root(root, true);

verifyEqual(testCase, actualRoot, root);
verifyTrue(testCase, is_folder_empty(root));
end

function testPreservesResultsWhenCleaningIsDisabled(testCase)
root = tempname;
mkdir(root);
cleanup = onCleanup(@() remove_if_present(root)); %#ok<NASGU>

oldResult = fullfile(root, 'keep.txt');
write_text(oldResult, 'keep');

actualRoot = hydropol2d_prepare_output_root(root, false);

verifyEqual(testCase, actualRoot, root);
verifyTrue(testCase, isfile(oldResult));
end

function write_text(path, value)
fileId = fopen(path, 'w');
assert(fileId >= 0, 'Could not create test file: %s', path);
cleanup = onCleanup(@() fclose(fileId)); %#ok<NASGU>
fwrite(fileId, value);
end

function tf = is_folder_empty(folder)
entries = dir(folder);
tf = all(ismember({entries.name}, {'.', '..'}));
end

function remove_if_present(folder)
if isfolder(folder)
    fileattrib(folder, '+w', 'u', 's');
    rmdir(folder, 's');
end
end
