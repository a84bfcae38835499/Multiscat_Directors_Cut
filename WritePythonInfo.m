function [] = WritePythonInfo(fileprefix,a1,a2,area,b1,b2,Nsuper,Theta,Nensemble,inputEntropy,defectDensity,Ndefect)
  writematrix([],fileprefix+'.info_for_vivian_python_nice_plotting_hexagon_script',FileType='text')
  a1str = [char(num2str(a1))];
  a2str = [char(num2str(a2))];
  b1str = [char(num2str(b1(1:2)))];
  b2str = [char(num2str(b2(1:2)))];
  areastr = [char(num2str(area))];
  nsupstr = [char(num2str(Nsuper))];
  thetastr = [char(num2str(Theta))];
  nens = [char(num2str(Nensemble))];
  denstr = [char(num2str(defectDensity))];
  entropstr = [char(num2str(inputEntropy,'%0.6f'))];
  defstr = [char(num2str(Ndefect))];
  S = fileread('latticeVects.info_for_vivian_python_nice_plotting_hexagon_script');
  realStr = ['Real space vectors:',newline,'a1 = ',a1str, newline, 'a2 = ',a2str,newline,'Unit cell area = ',areastr,newline,'Nsuper = ',nsupstr];
  recpStr = ['Reciprocal vectors:',newline,'b1 = ',b1str, newline, 'b2 = ', b2str];
  defectStr = ['Defect data:',newline,'Theta = ', thetastr,newline,'Ensemble size = ', ...
    nens,newline,'Positional entropy = ',entropstr,newline, ...
    'Defect density in cm^-2 = ',denstr,newline,'Total number of defects = ',defstr];
  S = [realStr,newline,recpStr,newline,defectStr,S];
  FID = fopen(fileprefix + '.info_for_vivian_python_nice_plotting_hexagon_script', 'w');
  if FID == -1, error('Cannot open file %s', FileName); end
  fwrite(FID, S, 'char');
  fclose(FID);
end