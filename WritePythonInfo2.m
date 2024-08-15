function [] = WritePythonInfo2(fileprefix,a1,a2,b1,b2,Nsuper)
  writematrix([],fileprefix+'.IFVPNPHS2',FileType='text')
  a1str = [char(num2str(a1))];
  a2str = [char(num2str(a2))];
  b1str = [char(num2str(b1(1:2)))];
  b2str = [char(num2str(b2(1:2)))];
  nsupstr = [char(num2str(Nsuper))];
  realStr = ['Real space vectors:',newline,'a1 = ',a1str, newline, 'a2 = ',a2str,newline,'Nsuper = ',nsupstr,newline];
  recpStr = ['Reciprocal vectors:',newline,'b1 = ',b1str, newline, 'b2 = ', b2str];
  S = [realStr,newline,recpStr];
  FID = fopen(fileprefix + '.IFVPNPHS2', 'w');
  if FID == -1, error('Cannot open file %s', FileName); end
  fwrite(FID, S, 'char');
  fclose(FID);
end