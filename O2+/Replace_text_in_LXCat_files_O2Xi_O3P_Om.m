% Replacing text in LXCat files for compatibility with LoKI-B
% 17.02.2021
fileID= fopen('Laporta_e_O2Xi__O3P_Om.txt','r');
fileOut= fopen('Cross_sections_out.txt','w');
% line = fgetl(fileID);
while true
    line = fgetl(fileID);
    if ~ischar(line) %eof
        break;
    end
    if length(line)>5
        switch line(1:5)
            case 'O2(v='
                num_lvl=line(6);
                if line(7)~=')'
                    num_lvl(2)=line(7);
                end
                line=['O2-elec -> O- + O (O2(X,v=' num_lvl ...
                                                    ') -> O- + O(3P))'];
            case 'SPECI'
                line='SPECIES: e / O2-elec';
            case 'PROCE'
                line=['PROCESS: E + O2-elec -> O- + O (O2(X,v=' num_lvl...
                    ') -> O- + O(3P)), Attachment'];
            case 'PARAM'
                fprintf(fileOut, [line '\n']);
                line=['COMMENT: [e + O2(X,v=' num_lvl ...
                    ') -> O(-,gnd) + O(3P), Attachment] Laporta 2015'];
        end
    end
    fprintf(fileOut, [line '\n']);
end
fclose all;