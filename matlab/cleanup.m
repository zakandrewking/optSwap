function cleanup()
   global run status 
   tw = setupTwitty;
   if ~isempty(tw)
       string = sprintf('@zakandrewking %s: %s', run, status);
       disp(string);
       string = string(1:min(140,length(string)));
       tw.updateStatus(string);
   else
       warning('Could not setup Twitty');
   end
end