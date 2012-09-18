function cleanup()
   global run status 
   tw = setupTwitty;
   string = sprintf('@zakandrewking %s: %s', run, status);
   disp(string);
   string = string(1:min(140,length(string)));
   tw.updateStatus(string);
end