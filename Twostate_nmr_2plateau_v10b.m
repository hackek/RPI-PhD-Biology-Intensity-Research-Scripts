%% THIS SCRIPT Fits equilibrium 2state pressure data when x(1) = 1 it is fraction native
clear all;
close all;
clc;

%prompt = 'Enter Temperature in Kelvin: ';
%Temp = input(prompt);
Temp = 298;
RGC = 0.00197;
RGP = 83.1451;
%prompt='Enter number of curves for entire experiment:';
%Num = input(prompt);
%prompt='Enter number of pressure points for each experiment:';
%Pres = input(prompt);

%Early attempts to read data
%T = readtable(filename,'ReadVariablenames',true);
%resname = T.Properties.VariableNames;
%resname2 = resname(1,2:(Num+1));

%This section opens and reads the data file.
%Data file must be an excel file with the first row empty. The second row
%must have the words Pressure in the first column and the residue numbers in
%the other columns. These numbers must be justified left with no spaces.
%Then rows 3-Pres are the values of the pressures (in column 1) and the
%reisdue peak intensities (or volumes as the case may be) for each pressure
%in the other columns.
%This section also figures out how many pressures and how many curves
%(residues) were measured and then ceastes files for outputting the results

%data_2023-09-19
prompt = 'Enter data filename (must be an .xlsx file, do not specify extension): ';
fid = input(prompt,'s');
data_filename = [fid, '.xlsx'];
[num,text,raw] = xlsread(data_filename);
%disp(num); %displays startng at row 2 all row/column data
%disp(text); %displays "{}(0x0)", which isn't used in the rest of the program
%disp(raw); %not a variable used in the program
dim = size(num);
%disp(dim); %starting at row2, displays inclusively # of rows and then # of columns
Pres = dim(1) - 1; %Gets the total rows minus the header row
Num = dim(2) - 1; %Gets the total columns minus the header column
resname2 = num(1,2:(Num+1));
%disp(resname2);%Gets the column headers in row 2

resfilename = [fid, 'r.xlsx']; % r file
fitparamfilename = [fid, 'p.xlsx']; %p file

%This section pre-allocates variables
IntF = zeros(Num);
sigma_IntF = zeros(Num);
IntU = zeros(Num);
sigma_IntU = zeros(Num);
DeltaG = zeros(Num);
sigma_DG = zeros(Num);
DeltaV = zeros(Num);
sigma_DV = zeros(Num);
fracintfit = zeros(Pres,Num);
pressureval= zeros(Pres,Num);
Ipcalcres = zeros(Num);
Ipcalcresn = zeros(Num);
Ipcalc = zeros(Pres);

% Because of memory issues with Matlab, I had to divide the analyses into
% chunks of 25 residues otherwise Matlab crashes. Num is the Number of
% measured residues, which is obtained from the number of columns of the
% datafile minus one (to subtract the column with the xvalues of pressure)
% Pres is the number of pressures also obtained from the number of rows in
% the data  file (minus % one row for the headers)

%This section figures out how many partitions there are %Numpart, Numdiv, Numleft
if Num < 25
    Numpart = 1;
    Numdiv = 1;
else
    Numdiv = Num./25;
    Numpart = 1;
    while Numpart <= Numdiv
        Numpart = Numpart + 1;
    end
    Numpart = Numpart - 1;
    Numleft = Num - (Numpart.*25);
end

% This section gets the parameter guesses
prompt='Enter the initial Ku (0.01 is a good guess): ';
Kuo=input(prompt);
prompt='Enter the intial deltaVu (50 is a good guess): ';
DVu=input(prompt);
DVu = DVu/(RGP.*Temp);
prompt='Set final plateau to 0 (Y or N)?: ';
str = input(prompt,'s'); %finalPlateau
if str == 'Y'
    IntHP = 0;
end
prompt = 'Fix initial plateau (Y or N)?: ';
str2 = input(prompt, 's'); %initialPlateau
    if str2 == 'Y'
        prompt='Fix all initial plateau values at 105% of average of first and second values (Y or N)?: ';
        str3 = input(prompt, 's');
    end
Part = 1;
col = 2;
count = 1;

%This section is used if there are fewer than 25 curves in the dataset
% Calculates and plots each point for every entry then assigns to r file
if Num < 25
    while count <= Num % For loop parsing over every column making these calculations
        colreal = ((Part - 1).* 25) + col;
        res = colreal - 1;

        xdata = num(2:Pres+1,1); % All row entries of column 1, assignmentValues
        ydata = num(2:Pres+1,colreal); % All row entries of column 2, curr values

        Int0 = (ydata(1) + ydata(2))./2;
        %Initial guess for low pressure plateau values is the intenity of
        %that residue in the first pressure point

        if str == 'N'
            IntHP= ydata(Pres);
        end

%If the fit will not be forced to a zero value at high pressure
%(str == N) then the initial guess for the high pressure plateau is
%the intensity value for that residue at the last pressure point

%Fit the data for residue (res)

        if str2 == 'Y'
            if str3 == 'N'
                f = figure; % defined but never used
                scatter(xdata,ydata);
                hold on;
                prompt='Enter initial plateau value to fix:   ';
                I0F = input(prompt);
                Int0 = I0F;
            else
                Int0 = Int0*1.05;
                I0F = Int0;
            end
        else
            I0F = 0;
        end

        if I0F == 0;
            if str == 'N'
               p0 = [Int0,IntHP,Kuo,DVu]; %function initial values = IntF, IntU, DG, DV
               fun = @(p,xdata) (p(1) + p(2).*p(3)*exp(p(4)*xdata))./(1 + (p(3)*exp(p(4)*xdata))); %matlab anonymous function
               format long; % float with 15 decimals after it
               %P, resid, and J are used, whereas resnorm, output, and lambda are not, aka (1, 3, 7) are used
               %1=x array the size of p0, 3=value of objective function at solution (fun(xdata)-ydata), 7=Jacobian matrix jac(i,j) partial div of fun(i) with respect to x(j)
               %1=x0 optimized values from p0, 3=fun(x0,xdata)-ydata optimized values array, 7=jac(i,j) matrix
               [P,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(fun,p0,xdata,ydata);
               ci = nlparci(P,resid,'jacobian',J); %95% confidence intervals w/coefficients, residuals, Jacovian
               confm_IntF = ci(1,1);%lower bound ci
               confp_IntF = ci(1,2);%upper bound ci
               confm_IntU = ci(2,1);
               confp_IntU = ci(2,2);
               confm_DG = ci(3,2);
               confp_DG = ci(3,2);
               confm_DV = ci(4,1);
               confp_DV = ci(4,2);
            else
                p0 = [Int0,Kuo,DVu];
                fun = @(p,xdata) (p(1).*(1./(1 + (p(2)*exp(p(3)*xdata)))));
                format long;
                [P,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(fun,p0,xdata,ydata);
                ci = nlparci(P,resid,'jacobian',J);
                confm_IntF = ci(1,1);
                confp_IntF = ci(1,2);
                confm_DG = ci(2,2);
                confp_DG = ci(2,2);
                confm_DV = ci(3,1);
                confp_DV = ci(3,2);
            end
        else
            if str == 'N'
               p0 = [IntHP,Kuo,DVu];
               fun = @(p,xdata) (I0F + p(1).*p(2)*exp(p(3)*xdata))./(1 + (p(2)*exp(p(3)*xdata)));
               format long;
               [P,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(fun,p0,xdata,ydata);
               ci = nlparci(P,resid,'jacobian',J);
               confm_IntU = ci(1,1);%lower bound ci
               confp_IntU = ci(1,2);%upper bound ci
               confm_DG = ci(2,2);
               confp_DG = ci(2,2);
               confm_DV = ci(3,1);
               confp_DV = ci(3,2);
            else
                p0 = [Kuo,DVu];
                fun = @(p,xdata) (I0F.*(1./(1 + (p(1)*exp(p(2)*xdata)))));
                format long;
                [P,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(fun,p0,xdata,ydata);
                ci = nlparci(P,resid,'jacobian',J);
                confm_DG = ci(1,1);
                confp_DG = ci(1,2);
                confm_DV = ci(2,1);
                confp_DV = ci(2,2);
            end
        end
%Plot the fit and data
        f=figure; % defined but never used
        plot(xdata,fun(P,xdata),'color','red');
        hold on;

        scatter(xdata,ydata);
        hold on

        legend(strcat('Fit', num2str(resname2(colreal-1))));
        title('High pressure dentauration');
        fid = num2str(resname2(colreal-1));
        plotfilename = [fid, '.tiff'];
        print(fullfile(pwd, plotfilename), '-dtiff', '-r300');
        %print(plotfilename, '-dpng');

%calculate the y-values from the fit for each pressure
        p = 1;
        while p <= Pres
            if I0F == 0;
                if str == 'N'
                    Ipcalc(p) = (P(1) + P(2).*P(3)*exp(P(4)*xdata(p)))./(1 + (P(3)*exp(P(4)*xdata(p))));
                else
                    Ipcalc(p) = (P(1).*(1./(1 + (P(2)*exp(P(3)*xdata(p))))));
                end
            else
                if str == 'N'
                    Ipcalc(p) = (I0F + P(1).*P(2)*exp(P(3)*xdata(p)))./(1 + (P(2)*exp(P(3)*xdata(p))));
                else
                    Ipcalc(p) = (I0F.*(1./(1 + (P(1)*exp(P(2)*xdata(p))))));
                end
            end
            p = p + 1;
        end
        Ipcalc2 =Ipcalc(:,1);

%Calculate the normalized data and fit for each pressure
%either with not forcing to 0 (str == N) or forcing to 0 (else)
        if I0F == 0
            if str == 'N'
                Ipcalc2n = (Ipcalc2 - P(2))./(P(1)- P(2));
                ydatan = (ydata - P(2))./(P(1) - P(2));
                resid2 = transpose(resid);
            else
                Ipcalc2n = (Ipcalc2 - P(2))./P(1);
                ydatan = (ydata - P(2))./P(1);
                resid2 = transpose(resid);
            end
        else
            if str == 'N'
                Ipcalc2n = (Ipcalc2 - P(1))./(I0F - P(1));
                ydatan = (ydata - P(1))./(I0F - P(1));
                resid2 = transpose(resid);
            else
                Ipcalc2n = (Ipcalc2)./I0F;
                ydatan = (ydata)./I0F;
                resid2 = transpose(resid);
            end
        end

% make a table of the data, fit, normalized data and normalized fit
        pressure = xdata;
        data = ydata;
        fit = Ipcalc2;
        residuals = resid2;
        normdata = ydatan;
        normfit = Ipcalc2n;

%Write data and fit to file
        warning( 'off', 'MATLAB:xlswrite:AddSheet' );
        sheet = num2str(resname2(res));
        results = [pressure(:),data(:),fit(:),resid2(:),normdata(:),normfit(:)];
        xlRange = 'A1';
        xlswrite(resfilename,results,sheet,xlRange);

%Fill out the parameters for the current residue
        if I0F == 0
            if str == 'N'
                IntF(res) = P(1);
                sigma_IntF(res) = (confp_IntF - P(1))./2;
                IntU(res) = P(2);
                sigma_IntU(res) = (confp_IntU - P(2))./2;
                DeltaG(res) = (-1.*Temp.*RGC).*log(P(3));
                confp2_DG = ((-1.*Temp.*RGC).*log(confp_DG));
                sigma_DG(res) = (DeltaG(res) - confp2_DG)./2;
                DeltaV(res) = P(4).* RGP.*Temp;
                confm2_DV = (confm_DV).*RGP.*Temp;
                sigma_DV(res) = (DeltaV(res) - confm2_DV)./2;
            else
                IntF(res) = P(1);
                sigma_IntF(res) = confp_IntF - P(1);
                DeltaG(res) = (-1.*Temp.*RGC).*log(P(2));
                confp2_DG = ((-1.*Temp.*RGC).*log(confp_DG));
                sigma_DG(res) = (DeltaG(res) - confp2_DG)./2;
                DeltaV(res) = P(3).* RGP.*Temp;
                confm2_DV = (confm_DV).*RGP.*Temp;
                sigma_DV(res) = (DeltaV(res) - confm2_DV)./2;
            end
        else
            if str == 'N'
                IntU(res) = P(1);
                sigma_IntU(res) = (confp_IntU - P(1))./2;
                DeltaG(res) = (-1.*Temp.*RGC).*log(P(2));
                confp2_DG = ((-1.*Temp.*RGC).*log(confp_DG));
                sigma_DG(res) = (DeltaG(res) - confp2_DG)./2;
                DeltaV(res) = P(3).* RGP.*Temp;
                confm2_DV = (confm_DV).*RGP.*Temp;
                sigma_DV(res) = (DeltaV(res) - confm2_DV)./2;
            else
                DeltaG(res) = (-1.*Temp.*RGC).*log(P(1));
                confp2_DG = ((-1.*Temp.*RGC).*log(confp_DG));
                sigma_DG(res) = (DeltaG(res) - confp2_DG)./2;
                DeltaV(res) = P(2).* RGP.*Temp;
                confm2_DV = (confm_DV).*RGP.*Temp;
                sigma_DV(res) = (DeltaV(res) - confm2_DV)./2;
            end
        end

        col = col+1;
        count = count + 1;
    end
else

%Fit each curve to a two state model when there are more than 25 curves
    Part = 1;
    col = 2;
    while Part <= Numpart
        count = 1;
        while count <= 25
            colreal = ((Part - 1).* 25) + col;
            res = colreal - 1;
            xdata = num(2:Pres+1,1);
            ydata = num(2:Pres+1,colreal);

            Int0 = (ydata(1) + ydata(2))./2;
            if str == 'N'
                IntHP = ydata(Pres);
            end

        if str2 == 'Y'
            if str3 == 'N'
                f = figure;
                scatter(xdata,ydata);
                hold on;
                prompt='Enter initial plateau value to fix:   ';
                I0F = input(prompt);
                Int0 = I0F;
            else
                Int0 = Int0*1.05;
                I0F = Int0;
            end
        else
            I0F = 0;
        end

    if I0F == 0;
        if str == 'N'
           p0 = [Int0,IntHP,Kuo,DVu];
           fun = @(p,xdata) (p(1) + p(2).*p(3)*exp(p(4)*xdata))./(1 + (p(3)*exp(p(4)*xdata)));
           format long;
           [P,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(fun,p0,xdata,ydata);
           ci = nlparci(P,resid,'jacobian',J);
           confm_IntF = ci(1,1);
           confp_IntF = ci(1,2);
           confm_IntU = ci(2,1);
           confp_IntU = ci(2,2);
           confm_DG = ci(3,2);
           confp_DG = ci(3,2);
           confm_DV = ci(4,1);
           confp_DV = ci(4,2);
        else
            p0 = [Int0,Kuo,DVu];
            fun = @(p,xdata) (p(1).*(1./(1 + (p(2)*exp(p(3)*xdata)))));
            format long;
            [P,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(fun,p0,xdata,ydata);
            ci = nlparci(P,resid,'jacobian',J);
            confm_IntF = ci(1,1);
            confp_IntF = ci(1,2);
            confm_DG = ci(2,2);
            confp_DG = ci(2,2);
            confm_DV = ci(3,1);
            confp_DV = ci(3,2);
        end
    else
        if str == 'N'
           p0 = [IntHP,Kuo,DVu];
           fun = @(p,xdata) (I0F + p(1).*p(2)*exp(p(3)*xdata))./(1 + (p(2)*exp(p(3)*xdata)));
           format long;
           [P,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(fun,p0,xdata,ydata);
           ci = nlparci(P,resid,'jacobian',J);
           confm_IntU = ci(1,1);
           confp_IntU = ci(1,2);
           confm_DG = ci(2,2);
           confp_DG = ci(2,2);
           confm_DV = ci(3,1);
           confp_DV = ci(3,2);
        else
            p0 = [Kuo,DVu];
            fun = @(p,xdata) (I0F.*(1./(1 + (p(1)*exp(p(2)*xdata)))));
            format long;
            [P,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(fun,p0,xdata,ydata);
            ci = nlparci(P,resid,'jacobian',J);
            confm_DG = ci(1,1);
            confp_DG = ci(1,2);
            confm_DV = ci(2,1);
            confp_DV = ci(2,2);
        end
    end
%Plot the fit and data
        f=figure; % defined but never used
        plot(xdata,fun(P,xdata),'color','red');
        hold on;

        scatter(xdata,ydata);
        hold on

        legend(strcat('Fit', num2str(resname2(colreal-1))));
        title('High pressure dentauration');
        fid = num2str(resname2(colreal-1));
        plotfilename = [fid, '.tiff'];
        print(fullfile(pwd, plotfilename), '-dtiff', '-r300');

%calculate the y-values from the fit for each pressure
        p = 1;
        while p <= Pres
            if I0F == 0;
                if str == 'N'
                    Ipcalc(p) = (P(1) + P(2).*P(3)*exp(P(4)*xdata(p)))./(1 + (P(3)*exp(P(4)*xdata(p))));
                else
                    Ipcalc(p) = (P(1).*(1./(1 + (P(2)*exp(P(3)*xdata(p))))));
                end
            else
                if str == 'N'
                    Ipcalc(p) = (I0F + P(1).*P(2)*exp(P(3)*xdata(p)))./(1 + (P(2)*exp(P(3)*xdata(p))));
                else
                    Ipcalc(p) = (I0F.*(1./(1 + (P(1)*exp(P(2)*xdata(p))))));
                end
            end
            p = p + 1;
        end
        Ipcalc2 =Ipcalc(:,1);

%Calculate the normalized data and fit for each pressure
%either with not forcing to 0 (str == N) or forcing to 0 (else)
        if I0F == 0
            if str == 'N'
                Ipcalc2n = (Ipcalc2 - P(2))./(P(1)- P(2));
                ydatan = (ydata - P(2))./(P(1) - P(2));
                resid2 = transpose(resid);
            else
                Ipcalc2n = (Ipcalc2 - P(2))./P(1);
                ydatan = (ydata - P(2))./P(1);
                resid2 = transpose(resid);
            end
        else
            if str == 'N'
                Ipcalc2n = (Ipcalc2 - P(1))./(I0F - P(1));
                ydatan = (ydata - P(1))./(I0F - P(1));
                resid2 = transpose(resid);
            else
                Ipcalc2n = (Ipcalc2)./I0F;
                ydatan = (ydata)./I0F;
                resid2 = transpose(resid);
            end
        end

% make a table of the data, fit, normalized data and normalized fit
        pressure = xdata;
        data = ydata;
        fit = Ipcalc2;
        residuals = resid2;
        normdata = ydatan;
        normfit = Ipcalc2n;

%Write data and fit to file
        warning( 'off', 'MATLAB:xlswrite:AddSheet' );
        sheet = num2str(resname2(res));
        results = [pressure(:),data(:),fit(:),resid2(:),normdata(:),normfit(:)];
        xlRange = 'A1';
        xlswrite(resfilename,results,sheet,xlRange);

%Fill out the parameters for the current residue
        if I0F == 0
            if str == 'N'
                IntF(res) = P(1);
                sigma_IntF(res) = (confp_IntF - P(1))./2;
                IntU(res) = P(2);
                sigma_IntU(res) = (confp_IntU - P(2))./2;
                DeltaG(res) = (-1.*Temp.*RGC).*log(P(3));
                confp2_DG = ((-1.*Temp.*RGC).*log(confp_DG));
                sigma_DG(res) = (DeltaG(res) - confp2_DG)./2;
                DeltaV(res) = P(4).* RGP.*Temp;
                confm2_DV = (confm_DV).*RGP.*Temp;
                sigma_DV(res) = (DeltaV(res) - confm2_DV)./2;
            else
                IntF(res) = P(1);
                sigma_IntF(res) = confp_IntF - P(1);
                DeltaG(res) = (-1.*Temp.*RGC).*log(P(2));
                confp2_DG = ((-1.*Temp.*RGC).*log(confp_DG));
                sigma_DG(res) = (DeltaG(res) - confp2_DG)./2;
                DeltaV(res) = P(3).* RGP.*Temp;
                confm2_DV = (confm_DV).*RGP.*Temp;
                sigma_DV(res) = (DeltaV(res) - confm2_DV)./2;
            end
        else
            if str == 'N'
                IntU(res) = P(1);
                sigma_IntU(res) = (confp_IntU - P(1))./2;
                DeltaG(res) = (-1.*Temp.*RGC).*log(P(2));
                confp2_DG = ((-1.*Temp.*RGC).*log(confp_DG));
                sigma_DG(res) = (DeltaG(res) - confp2_DG)./2;
                DeltaV(res) = P(3).* RGP.*Temp;
                confm2_DV = (confm_DV).*RGP.*Temp;
                sigma_DV(res) = (DeltaV(res) - confm2_DV)./2;
            else
                DeltaG(res) = (-1.*Temp.*RGC).*log(P(1));
                confp2_DG = ((-1.*Temp.*RGC).*log(confp_DG));
                sigma_DG(res) = (DeltaG(res) - confp2_DG)./2;
                DeltaV(res) = P(2).* RGP.*Temp;
                confm2_DV = (confm_DV).*RGP.*Temp;
                sigma_DV(res) = (DeltaV(res) - confm2_DV)./2;
            end
        end

        col = col+1;
        count = count + 1;
        end
    col = 2;
    Part = Part + 1;

    end
    col = 2;
    count = 1;
%This part fits the curves remaining after the last partiton into 25
    while count <= Numleft
        colreal = ((Numpart).* 25) + (count + 1);
        res = colreal - 1;

        xdata = num(2:Pres+1,1);
        ydata = num(2:Pres+1,colreal);

        Int0 = (ydata(1) + ydata(2))./2;
        if str == 'N'
            IntHP = ydata(Pres);
        end

        if str2 == 'Y'
            if str3 == 'N'
                f = figure;
                scatter(xdata,ydata);
                hold on;
                prompt='Enter initial plateau value to fix:   ';
                I0F = input(prompt);
                Int0 = I0F;
            else
                Int0 = Int0*1.05;
                I0F = Int0;
            end
        else
            I0F = 0;
        end

    if I0F == 0;
        if str == 'N'
           p0 = [Int0,IntHP,Kuo,DVu];
           fun = @(p,xdata) (p(1) + p(2).*p(3)*exp(p(4)*xdata))./(1 + (p(3)*exp(p(4)*xdata)));
           format long;
           [P,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(fun,p0,xdata,ydata);
           ci = nlparci(P,resid,'jacobian',J);
           confm_IntF = ci(1,1);
           confp_IntF = ci(1,2);
           confm_IntU = ci(2,1);
           confp_IntU = ci(2,2);
           confm_DG = ci(3,2);
           confp_DG = ci(3,2);
           confm_DV = ci(4,1);
           confp_DV = ci(4,2);
        else
            p0 = [Int0,Kuo,DVu];
            fun = @(p,xdata) (p(1).*(1./(1 + (p(2)*exp(p(3)*xdata)))));
            format long;
            [P,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(fun,p0,xdata,ydata);
            ci = nlparci(P,resid,'jacobian',J);
            confm_IntF = ci(1,1);
            confp_IntF = ci(1,2);
            confm_DG = ci(2,2);
            confp_DG = ci(2,2);
            confm_DV = ci(3,1);
            confp_DV = ci(3,2);
        end
    else
        if str == 'N'
           p0 = [IntHP,Kuo,DVu];
           fun = @(p,xdata) (I0F + p(1).*p(2)*exp(p(3)*xdata))./(1 + (p(2)*exp(p(3)*xdata)));
           format long;
           [P,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(fun,p0,xdata,ydata);
           ci = nlparci(P,resid,'jacobian',J);
           confm_IntU = ci(1,1);
           confp_IntU = ci(1,2);
           confm_DG = ci(2,2);
           confp_DG = ci(2,2);
           confm_DV = ci(3,1);
           confp_DV = ci(3,2);
        else
            p0 = [Kuo,DVu];
            fun = @(p,xdata) (I0F.*(1./(1 + (p(1)*exp(p(2)*xdata)))));
            format long;
            [P,resnorm,resid,exitflag,output,lambda,J] = lsqcurvefit(fun,p0,xdata,ydata);
            ci = nlparci(P,resid,'jacobian',J);
            confm_DG = ci(1,1);
            confp_DG = ci(1,2);
            confm_DV = ci(2,1);
            confp_DV = ci(2,2);
        end
    end
%Plot the fit and data
        f=figure; % defined but never used
        plot(xdata,fun(P,xdata),'color','red');
        hold on;

        scatter(xdata,ydata);
        hold on

        legend(strcat('Fit', num2str(resname2(colreal-1))));
        title('High pressure dentauration');
        fid = num2str(resname2(colreal-1));
        plotfilename = [fid, '.tiff'];
        print(fullfile(pwd, plotfilename), '-dtiff', '-r300');

%calculate the y-values from the fit for each pressure
        p = 1;
        while p <= Pres
            if I0F == 0;
                if str == 'N'
                    Ipcalc(p) = (P(1) + P(2).*P(3)*exp(P(4)*xdata(p)))./(1 + (P(3)*exp(P(4)*xdata(p))));
                else
                    Ipcalc(p) = (P(1).*(1./(1 + (P(2)*exp(P(3)*xdata(p))))));
                end
            else
                if str == 'N'
                    Ipcalc(p) = (I0F + P(1).*P(2)*exp(P(3)*xdata(p)))./(1 + (P(2)*exp(P(3)*xdata(p))));
                else
                    Ipcalc(p) = (I0F.*(1./(1 + (P(1)*exp(P(2)*xdata(p))))));
                end
            end
            p = p + 1;
        end
        Ipcalc2 =Ipcalc(:,1);

%Calculate the normalized data and fit for each pressure
%either with not forcing to 0 (str == N) or forcing to 0 (else)
        if I0F == 0
            if str == 'N'
                Ipcalc2n = (Ipcalc2 - P(2))./(P(1)- P(2));
                ydatan = (ydata - P(2))./(P(1) - P(2));
                resid2 = transpose(resid);
            else
                Ipcalc2n = (Ipcalc2 - P(2))./P(1);
                ydatan = (ydata - P(2))./P(1);
                resid2 = transpose(resid);
            end
        else
            if str == 'N'
                Ipcalc2n = (Ipcalc2 - P(1))./(I0F - P(1));
                ydatan = (ydata - P(1))./(I0F - P(1));
                resid2 = transpose(resid);
            else
                Ipcalc2n = (Ipcalc2)./I0F;
                ydatan = (ydata)./I0F;
                resid2 = transpose(resid);
            end
        end

% make a table of the data, fit, normalized data and normalized fit
        pressure = xdata;
        data = ydata;
        fit = Ipcalc2;
        residuals = resid2;
        normdata = ydatan;
        normfit = Ipcalc2n;

%Write data and fit to file
        warning( 'off', 'MATLAB:xlswrite:AddSheet' );
        sheet = num2str(resname2(res));
        results = [pressure(:),data(:),fit(:),resid2(:),normdata(:),normfit(:)];
        xlRange = 'A1';
        xlswrite(resfilename,results,sheet,xlRange);

%Fill out the parameters for the current residue
        if I0F == 0
            if str == 'N'
                IntF(res) = P(1);
                sigma_IntF(res) = (confp_IntF - P(1))./2;
                IntU(res) = P(2);
                sigma_IntU(res) = (confp_IntU - P(2))./2;
                DeltaG(res) = (-1.*Temp.*RGC).*log(P(3));
                confp2_DG = ((-1.*Temp.*RGC).*log(confp_DG));
                sigma_DG(res) = (DeltaG(res) - confp2_DG)./2;
                DeltaV(res) = P(4).* RGP.*Temp;
                confm2_DV = (confm_DV).*RGP.*Temp;
                sigma_DV(res) = (DeltaV(res) - confm2_DV)./2;
            else
                IntF(res) = P(1);
                sigma_IntF(res) = confp_IntF - P(1);
                DeltaG(res) = (-1.*Temp.*RGC).*log(P(2));
                confp2_DG = ((-1.*Temp.*RGC).*log(confp_DG));
                sigma_DG(res) = (DeltaG(res) - confp2_DG)./2;
                DeltaV(res) = P(3).* RGP.*Temp;
                confm2_DV = (confm_DV).*RGP.*Temp;
                sigma_DV(res) = (DeltaV(res) - confm2_DV)./2;
            end
        else
            if str == 'N'
                IntU(res) = P(1);
                sigma_IntU(res) = (confp_IntU - P(1))./2;
                DeltaG(res) = (-1.*Temp.*RGC).*log(P(2));
                confp2_DG = ((-1.*Temp.*RGC).*log(confp_DG));
                sigma_DG(res) = (DeltaG(res) - confp2_DG)./2;
                DeltaV(res) = P(3).* RGP.*Temp;
                confm2_DV = (confm_DV).*RGP.*Temp;
                sigma_DV(res) = (DeltaV(res) - confm2_DV)./2;
            else
                DeltaG(res) = (-1.*Temp.*RGC).*log(P(1));
                confp2_DG = ((-1.*Temp.*RGC).*log(confp_DG));
                sigma_DG(res) = (DeltaG(res) - confp2_DG)./2;
                DeltaV(res) = P(2).* RGP.*Temp;
                confm2_DV = (confm_DV).*RGP.*Temp;
                sigma_DV(res) = (DeltaV(res) - confm2_DV)./2;
            end
        end

        col = col+1;
        count = count + 1;
    end
end
res = 1;
while res <= Num % Makes sure all values if NAN are instead 0
    %not needed for Python
    residue(res) = res;
    if I0F == 0
       if str == 'N'
            IntF(isnan(IntF))=0;
            sigma_IntF(isnan(sigma_IntU)) = 0;
            IntU(isnan(IntU))=0;
            sigma_IntU(isnan(sigma_IntU)) = 0;
            DeltaG(isnan(DeltaG)) = 0;
            sigma_DG(isnan(sigma_DG)) = 0;
            DeltaV(isnan(DeltaV)) = 0;
            sigma_DV(isnan(sigma_DV)) = 0;
       else
            IntF(isnan(IntF))=0;
            sigma_IntF(isnan(sigma_IntF)) = 0;
            DeltaG(isnan(DeltaG)) = 0;
            sigma_DG(isnan(sigma_DG)) = 0;
            DeltaV(isnan(DeltaV)) = 0;
            sigma_DV(isnan(sigma_DV)) = 0;
       end
    else
       if str == 'N'
            IntU(isnan(IntU))=0;
            sigma_IntU(isnan(sigma_IntU)) = 0;
            DeltaG(isnan(DeltaG)) = 0;
            sigma_DG(isnan(sigma_DG)) = 0;
            DeltaV(isnan(DeltaV)) = 0;
            sigma_DV(isnan(sigma_DV)) = 0;
        else
            DeltaG(isnan(DeltaG)) = 0;
            sigma_DG(isnan(sigma_DG)) = 0;
            DeltaV(isnan(DeltaV)) = 0;
            sigma_DV(isnan(sigma_DV)) = 0;
        end
    end
    res = res + 1;
end

% Here's where the final data is assigned then writen to the p file
    resname3 = transpose(resname2); % Column headers flipped to be vertical
    res = res - 1; %Seems to have no effect as res is not referrenced here onwards
    residuenum = transpose(residue); %Seems to have no effect as res is not referrenced here onwards
    if I0F == 0
        if str == 'N'
            Int_Fold = IntF(:,1);
            sigma_IntF = sigma_IntF(:,1);
            Int_Unf = IntU(:,1);
            sigma_IntU = sigma_IntU(:,1);
            DeltaGu = DeltaG(:,1);
            sigma_DG = sigma_DG(:,1);
            DeltaVu = DeltaV(:,1);
            sigma_DV = sigma_DV(:,1);
            T2 = table(resname3,Int_Fold,sigma_IntF,Int_Unf,sigma_IntU,DeltaGu,sigma_DG,DeltaVu,sigma_DV);
            writetable(T2,fitparamfilename);
        else
            Int_Fold = IntF(:,1);
            sigma_IntF = sigma_IntF(:,1);
            DeltaGu = DeltaG(:,1);
            sigma_DG = sigma_DG(:,1);
            DeltaVu = DeltaV(:,1);
            sigma_DV = sigma_DV(:,1);
            T2 = table(resname3,Int_Fold,sigma_IntF,DeltaGu,sigma_DG,DeltaVu,sigma_DV);
            writetable(T2,fitparamfilename);
        end
    else
        if str == 'N'
            Int_Unf = IntU(:,1);
            sigma_IntU = sigma_IntU(:,1);
            DeltaGu = DeltaG(:,1);
            sigma_DG = sigma_DG(:,1);
            DeltaVu = DeltaV(:,1);
            sigma_DV = sigma_DV(:,1);
            T2 = table(resname3,Int_Fold,Int_Unf,sigma_IntU,DeltaGu,sigma_DG,DeltaVu,sigma_DV);
            writetable(T2,fitparamfilename);
        else
            DeltaGu = DeltaG(:,1);
            sigma_DG = sigma_DG(:,1);
            DeltaVu = DeltaV(:,1);
            sigma_DV = sigma_DV(:,1);
            T2 = table(resname3,DeltaGu,sigma_DG,DeltaVu,sigma_DV);
            writetable(T2,fitparamfilename);
        end
    end
    clear all:
    close all;
