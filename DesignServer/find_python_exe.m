%Find Python executable 
function pyexe = find_python_exe()
    pe = pyenv;
    if pe.Version ~= ""
            pyexe = pe.Executable;
            stderr = 0;
    else
        if ismac
            [stderr, pyexe] = system('which python3');
        elseif isunix
            [stderr, pyexe] = system('which python');
        elseif ispc
            [stderr, pyexe] = system('where python');
        else
            stderr = 1;
        end        
    end
    if stderr == 0
        pyexe = strip(pyexe);
    else
        pyexe = "";
        disp('Unable to find Python executable')
    end
end