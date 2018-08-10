function main
clearvars
clc
%Fucntion main calls the function _____
%Function (Directory of forward reads, Directory of reverse reads, Name of
%output file, parameter file).

%Local injection samples
nexterasort('input\RUN12\S22_S22_L001_R1_001.fastq\S22_S22_L001_R1_001.fastq','input\RUN12\S22_S22_L001_R2_001.fastq\S22_S22_L001_R2_001.fastq','output\S22GRNA2','pargRNA2.dat');
w=waitforbuttonpress();
end

function nexterasort(fwdread,revread,outfile1,parfile)
forwardreads = fastqread(fwdread);  %Read in forward reads
reversereads = fastqread(revread);  %Read in reverse reads

%Read in experiment specific parameters-------------------------------------------------------------------
fid = fopen(parfile);                   %Open parameter file
Intact = fgetl(fid);                    %Line 1 - Amplicon sequence if no editing takes place
Deletion = fgetl(fid);                  %Line 2 - Amplicon sequence if a targeted deletion occurs
Inversion = fgetl(fid);                 %Line 3 - Amplicon sequence if a targeted inversion occurs
DSBbp = str2num(fgetl(fid));            %Line 4 - DSB location in bp
WINDOWbp = str2num(fgetl(fid));         %Line 5 - BP window to check for indels on each side of DSB
MINSCORE = str2num(fgetl(fid));         %Line 6 - Minimum score for aligning 
SHORTMINSCORE = str2num(fgetl(fid));    %Line 7 - Minimum score for aligning short portion
AAVSCORE = str2num(fgetl(fid));     	%Line 8 - Minimum alignment score for the AAV genome
dysinv =  fgetl(fid);                   %Line 9 - Large portion of target gene with inversion for unexpected modifications (e.g. chewback)
dysrev = fgetl(fid);                    %Line 10 -  Large portion of target gene for unexpected modifications (e.g. chewback)
minsize = str2num(fgetl(fid));          %Line 11 - minimum amplicon size. Reject shorter amplicons
AAVGENOME1 = fgetl(fid);                %Line 12 - AAV Genome 1
AAVGENOME2 = fgetl(fid);                %Line 13 - AAV genome 2
DYSGapOpenValue = str2num(fgetl(fid));  %Line 14 - Value for gap opening when aligning to the dystrophin gene
AAVGAPOPENVALUE = str2num(fgetl(fid));  %Line 15 - Value for gap opening when aligning to the AAV genome
DYSEXTENDGAP = str2num(fgetl(fid));     %Line 16 - Value for gap extending for AAV genome
fclose(fid);
%END FILE READ----------------------------------------------------------

%Define remaining variables-----------------------------------------------
llim = 1;  %Start bp of amplicon is 1
winleft = DSBbp-WINDOWbp; %Window around DSB
winright = DSBbp+WINDOWbp; %Window around DSB
newStr = extractBetween(Intact,winleft+1,winright);
intexact = newStr{1};
newStr = extractBetween(Deletion,winleft+1,winright);
delexact = newStr{1};
newStr = extractBetween(Inversion,winleft+1,winright);
invexact = newStr{1};
newStr = extractBetween(Intact,1,DSBbp);
LeftAlign = newStr{1}; %Get DNA sequence from the start of the primer to the location of the DSB
AAVGENOME1R = seqrcomplement(AAVGENOME1);
AAVGENOME2R = seqrcomplement(AAVGENOME1);
%--------------------------------------------------------------------------

%initialize Variables------------------------------------------------------
delcount = 0;       %Number of deletions reads
invcount = 0;       %Number of inversions reads
intcount = 0;       %Number of intact reads
aavcount = 0;       %Number of AAV integration events
mysterycount = 0;   %Reads that did not align
i = 0;              %Counter for iterating through all reads and placing in matlab structure
j = 0;              %Counter for iterating through all reads and aligning
intactsubcount = 0; %Count the intact reads that have a substituion within the window
intactNHEJcount=0;  %count the intact reads that have an insertion/deletion/substitution
delindelcount = 0;  %Count the deletion reads that have an insertion or deletion
invindelcount = 0;  %Count the inverion reads that have an insertion or deletion
intactindelcount=0; %Count indels in intact reads

intactmatrix = cell(4000,1);    %output table for intact reads
deletionmatrix = cell(1000,1);  %output table for deletion reads
inversionmatrix = cell(1000,1); %Output table for inverted reads
aavintmatrix = cell(1000,1);    %Output table for aav integration reads
mysterymatrix = cell(4000,1);   %Output table for unalinged reads
indelmatrix = cell(4000,1);     %output table for indels
delindelmatrix = cell(4000,1);  %output table for deletion-indels
invindelmatrix = cell(4000,1);  %output table for inversion-indels


%Convert fastq to matlab table stuctures-----------------
structsize = size(forwardreads);
revstructsize = size(reversereads);
for i = 1:structsize(1,2)
    forwardreadS(i).Sequence = forwardreads(i).Sequence;
end
for rccount = 1:revstructsize(1,2)
    reversereadS(rccount).Sequence = seqrcomplement(reversereads(rccount).Sequence);
end
T = struct2table(forwardreadS);
T2 = struct2table(reversereadS);
%------------------------------------------------------

%Align each read using counter j
for j = 1:height(T)
    temp = T{j,1};  %Forward read
    S = char(temp);
    temp2 = T2{j,1}; %Reverse read
    S2 = char(temp2);
       
    if(cellfun('length',temp) > minsize) %Skip short amplicons
        [s0,Alignment]= swalign(S(llim:DSBbp), LeftAlign);
               
        if(s0 > SHORTMINSCORE)  %Remove reads that don't align to given locus
            
            [sf, af] = nwalign(dysrev,S(llim:length(S)), 'GapOpen', DYSGapOpenValue, 'Glocal', true,'ExtendGap', DYSEXTENDGAP);
            [sr, ar] = nwalign(dysrev,S2(llim:length(S2)), 'GapOpen', DYSGapOpenValue, 'Glocal', true,'ExtendGap', DYSEXTENDGAP);
            [si, ai] = nwalign(dysinv,S2(llim:length(S2)), 'GapOpen', DYSGapOpenValue, 'Glocal', true,'ExtendGap', DYSEXTENDGAP);
            Arr = strfind(af(3,:),'A');
            Trr = strfind(af(3,:),'T');
            Crr = strfind(af(3,:),'C');
            Grr = strfind(af(3,:),'G');
            tempstart = min([Arr,Trr,Crr,Grr]);
            tempend = max([Arr,Trr,Crr,Grr]);
            Arr = strfind(ar(3,:),'A');
            Trr = strfind(ar(3,:),'T');
            Crr = strfind(ar(3,:),'C');
            Grr = strfind(ar(3,:),'G');
            tempstart2 = min([Arr,Trr,Crr,Grr]);
            tempend2 = max([Arr,Trr,Crr,Grr]);
            Arr = strfind(ai(3,:),'A');
            Trr = strfind(ai(3,:),'T');
            Crr = strfind(ai(3,:),'C');
            Grr = strfind(ai(3,:),'G');
            tempstart3 = min([Arr,Trr,Crr,Grr]);
            tempend3 = max([Arr,Trr,Crr,Grr]);
            
            %align to reference and discard poor scores
			[s1, AlignInt] = swalign(S, Intact);
            [s2, Alignment] = swalign(S, Deletion);
		    [s3, Alignment] = swalign(S, Inversion);
          
           
           %CONDITION 1 - Read is an intact read-------------------------------------------------------------
            if (s1>s2)&&(s1>s3)&&(s1>MINSCORE)
                                              
                %PERFECT INTACT
                if(strfind(S,intexact)< 151)  
                     intcount=intcount+1;
                     s=strcat(S,':',S2,':intact:',num2str(tempstart),':',num2str(tempend),':',num2str(tempstart2),':',num2str(tempend2));
                     intactmatrix(intcount) = cellstr(s);
                %INDELS AND SUBS
                else
                    %CHECK FOR DASHES
                    DashFIND = strfind(AlignInt(1,:),'-');
                    DashFIND2 = strfind(AlignInt(3,:),'-');
                    DashFINDlength = length(DashFIND);
                    DashFINDlength2 = length(DashFIND2);
                    %SUBSITUTION
                    if(DashFINDlength==0)&&(DashFINDlength2==0)
                        intactNHEJcount = intactNHEJcount + 1;
                        intactsubcount = intactsubcount + 1;
                        s=strcat(S,':',S2,':SUB:',num2str(tempstart),':',num2str(tempend),':',num2str(tempstart2),':',num2str(tempend2));
                        indelmatrix(intactNHEJcount) = cellstr(s);
                   %INSERTION OR DELETION
                    else
                        %record indel in table and count
                        intactNHEJcount = intactNHEJcount + 1;
                        if ((min([DashFIND 151])<winright)&&(max([DashFIND 1])>winleft)) || ((min([DashFIND2 151])<winright)&&(max([DashFIND2 1])>winleft))
                            s=strcat(S,':',S2,':Arealindel:',num2str(tempstart),':',num2str(tempend),':',num2str(tempstart2),':',num2str(tempend2),':',num2str(DashFINDlength),':',num2str(DashFINDlength2),':',num2str(min(DashFIND)),':',num2str(max(DashFIND)),':',num2str(min(DashFIND2)),':',num2str(max(DashFIND2)));
                            intactindelcount = intactindelcount + 1;
                        else %out of window
                            s=strcat(S,':',S2,':outsideindel:',num2str(tempstart),':',num2str(tempend),':',num2str(tempstart2),':',num2str(tempend2),':',num2str(DashFINDlength),':',num2str(DashFINDlength2),':',num2str(min(DashFIND)),':',num2str(max(DashFIND)),':',num2str(min(DashFIND2)),':',num2str(max(DashFIND2)));
                        end    
                        indelmatrix(intactNHEJcount) = cellstr(s);
                    end
                end
            %END INTACT----------------------------------------------------
                    
            %CONDITION 2 - Read is a Deletion------------------------------
            elseif(s2>s3)&&(s2>MINSCORE)
                delcount=delcount+1;
                s=strcat(S,':',S2,':del:',num2str(tempstart),':',num2str(tempend),':',num2str(tempstart2),':',num2str(tempend2));
                deletionmatrix(delcount) = cellstr(s);
                if(strfind(S,delexact)< 151) 
                      %do nothing
                else
                    %record indel in table and count
                    delindelcount = delindelcount + 1;
                    delindelmatrix(delindelcount) = cellstr(S);
                end
            %END DELETION -------------------------------------------------
            
            %CONDITION 3 - Read is an Inversion----------------------------
            elseif(s3>MINSCORE)
                invcount=invcount+1;
                s=strcat(S,':',S2,':inv:',num2str(tempstart),':',num2str(tempend),':',num2str(tempstart3),':',num2str(tempend3));
                inversionmatrix(invcount) = cellstr(s);
                if(strfind(S,invexact)< 151) 
                      %do nothing
                else
                    %record indel in table and count
                    invindelcount = invindelcount + 1;
                    invindelmatrix(invindelcount) = cellstr(S);
                end
            %END INVERSION
            
            %include AAV integration etc 
            else
                [s4, a4] = nwalign(AAVGENOME1,S(DSBbp:length(S)), 'GapOpen', AAVGAPOPENVALUE,'Glocal', true);
                [s5, a5] = nwalign(AAVGENOME1R,S(DSBbp:length(S)), 'GapOpen', AAVGAPOPENVALUE,'Glocal', true);
                [s6, a6] = nwalign(AAVGENOME2,S(DSBbp:length(S)), 'GapOpen', AAVGAPOPENVALUE,'Glocal', true);
                [s7, a7] = nwalign(AAVGENOME2R,S(DSBbp:length(S)), 'GapOpen', AAVGAPOPENVALUE,'Glocal', true);

                if (s4>s5)&&(s4>s6)&&(s4>s7)&&(s4>AAVSCORE)
                    aavcount=aavcount+1;
                    Arr = strfind(a4(3,:),'A');
                    Trr = strfind(a4(3,:),'T');
                    Crr = strfind(a4(3,:),'C');
                    Grr = strfind(a4(3,:),'G');
                    temploc = min([Arr,Trr,Crr,Grr]);
                    templocMax = max([Arr,Trr,Crr,Grr]);
                    s=strcat(S,':',S2,':AAVGENOME1:',num2str(temploc),':',num2str(templocMax));
                    aavintmatrix(aavcount) = cellstr(s);

                elseif (s5>s6)&&(s5>s7)&&(s5>AAVSCORE)
                    aavcount=aavcount+1;
                    Arr = strfind(a5(3,:),'A');
                    Trr = strfind(a5(3,:),'T');
                    Crr = strfind(a5(3,:),'C');
                    Grr = strfind(a5(3,:),'G');
                    temploc = min([Arr,Trr,Crr,Grr]);
                    templocMax = max([Arr,Trr,Crr,Grr]);
                    s=strcat(S,':',S2,':AAVGENOME1R:',num2str(temploc),':',num2str(templocMax));
                    aavintmatrix(aavcount) = cellstr(s);
                   
                elseif (s6>s7)&&(s6>AAVSCORE)
                    aavcount=aavcount+1;
                    Arr = strfind(a6(3,:),'A');
                    Trr = strfind(a6(3,:),'T');
                    Crr = strfind(a6(3,:),'C');
                    Grr = strfind(a6(3,:),'G');
                    temploc = min([Arr,Trr,Crr,Grr]);
                    templocMax = max([Arr,Trr,Crr,Grr]);
                    s=strcat(S,':',S2,':AAVGENOME2:',num2str(temploc),':',num2str(templocMax));
                    aavintmatrix(aavcount) = cellstr(s);

                elseif(s7>AAVSCORE)
                    aavcount=aavcount+1;
                    Arr = strfind(a7(3,:),'A');
                    Trr = strfind(a7(3,:),'T');
                    Crr = strfind(a7(3,:),'C');
                    Grr = strfind(a7(3,:),'G');
                    temploc = min([Arr,Trr,Crr,Grr]);
                    templocMax = max([Arr,Trr,Crr,Grr]);
                    s=strcat(S,':',S2,':AAVGENOME2R:',num2str(temploc),':',num2str(templocMax));
                    aavintmatrix(aavcount) = cellstr(s);
                    
                else  %Try to align to the entire dystrophin gene
                    [catchscore, Alignment] = nwalign(dysrev,S(1:length(S)), 'GapOpen', DYSGapOpenValue,'Glocal', true,'ExtendGap', DYSEXTENDGAP);
                    if(catchscore >300)
                        mysterycount = mysterycount+1;
                        mysterymatrix(mysterycount) = cellstr(strcat(S,'DYSTROPHIN'));
                    else %Last - assign to mystery for manual analysis
                        mysterycount = mysterycount+1;
                        mysterymatrix(mysterycount) = cellstr(S);  
                    end
                 end
            end
        else %Do nothing if 5' end of amplicon doesn't align to target gene
        end
    else%Do nothing if amplicon too short
    end
end

%Output all numerical data and aligned reads into separate files
T = table(intcount,intactindelcount,delcount,invcount,aavcount,mysterycount)
%intcount
%intactNHEJcount
%intactsubcount
%delcount
%delindelcount
%invcount
%invindelcount
%aavcount
%mysterycount

fileID = fopen('output\allquantitative.dat','a'); %append data to outputfile
fprintf(fileID,'%6d\t %6d\t %6d\t %6d\t %6d\t %6d\t %6d\t %6d\t %6d\t %6d\n',[intcount intactindelcount intactNHEJcount intactsubcount delcount delindelcount invcount invindelcount aavcount mysterycount]);

Tabdata1 = cell2table(mysterymatrix);
writetable(Tabdata1,strcat(outfile1,'_mys.dat'))
Tabdata2 = cell2table(aavintmatrix);
writetable(Tabdata2,strcat(outfile1,'_aav.dat'))
Tabdata3 = cell2table(deletionmatrix);
writetable(Tabdata3,strcat(outfile1,'_del.dat'))
Tabdata4 = cell2table(indelmatrix);
writetable(Tabdata4,strcat(outfile1,'_indel.dat'))
Tabdata5 = cell2table(inversionmatrix);
writetable(Tabdata5,strcat(outfile1,'_inv.dat'))
Tabdata6 = cell2table(intactmatrix);
writetable(Tabdata6,strcat(outfile1,'_intact.dat'))
end
