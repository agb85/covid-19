function m = covidMobModel(tVec, bMatrix, M0, modelType)
% rI0 = 1;
    % Sabina/DMB mobility model (March 2021)
    % Age-structured COVID-19 model with mobility
    
    % Inputs: 
    % tVec = 24:14:(day of last estimated beta), 
    % bMatrix = a matrix with Nr rows and length(tVec) columns where the 
    %   (i, j)th entry is an initial guess for the jth beta in region i
    % M0 is the mobility scaling
    % modelType = 1 for metro, and 2 for lpha
    
    %% Set up the model with population and contact data
    I0 = 1; % Initial infection
    evtBd = 1; 
    tVec = [10, tVec]; % day zero is 12/31/2019 (this is an artifact of how I was incorporating mobility previously)
    tVec = tVec - 1;
    
    load('popInfo.mat', 'lpha1', 'lpha2', 'lphaByAge', 'metro', 'metroByAge')
    load('contactMatrices.mat', 'lphaContact', 'metroContact')
    load('vaccMatrices.mat', 'lphaVacMat', 'metVacMat')

    if (modelType == 1)
        
        % Jefferson includes Gilpin and Clear Creek
        regions = metro;
        popByAge = metroByAge;
        contactMat = metroContact; 
        vacMat = metVacMat;
          
    elseif (modelType == 2)
        
        regions = lpha2;
        popByAge = lphaByAge;
        contactMat = lphaContact;
        vacMat = lphaVacMat;
        
    end
    
    while size(contactMat, 3) < (length(tVec) + 1)
        % if simulation runs past available mobility data, continue using last available 
        contactMat(:, :, (size(contactMat, 3) + 1)) = contactMat(:, :, size(contactMat, 3));
    end
    regVec = 1:numel(regions);
    disp(regions)
    
    ageCat = {'0to19','20to39','40to64','65toinf'};
    Nr = length(regVec);
    Na = length(ageCat); % 4 age compartments
    
    m = sbiomodel('covidMobModel');   

    %% Parameters
    % Global constants
    p(1) = addparameter(m, 'beta', 0.4793); % transmission rate
    p(2) = addparameter(m, 'lambda', 1.395); % difference in infectiousness: symptomatic/asymptomatic
    p(3) = addparameter(m, 'alpha', 4); % incubation period
    p(4) = addparameter(m, 'gamma', 1/9); % recovery rate
    p(5) = addparameter(m, 'I0', I0); % initial number of infected individuals
    p(6) = addparameter(m, 'one', 1);
    p(7) = addparameter(m, 't0', 1);
    p(8) = addparameter(m, 'M', M0); % set(p(8), 'ConstantValue', false); % mobility scaling parameter
    
    % Parameter value switch times
    p(10) = addparameter(m, 'ioBegin', 24 - 1); % initial infection begins day 24 = 01/24/2020
    p(11) = addparameter(m, 'ioBd', 24 + evtBd - 1);
    p(12) = addparameter(m, 'vaxOn', 364 - 1); % vaccination begins December 15, move to V compartment after 14 days (Dec 29 = tv1 + 14)
    p(13) = addparameter(m, 'vaxS1', 382 - 1); % (Jan 16 = tv1 + 32)
    p(14) = addparameter(m, 'vaxS2', 412 - 1); % (Feb 15 = tv2 + 14)
    p(15) = addparameter(m, 'vaxS3', 430 - 1); % (March 5 = tv2 + 32)
    p(16) = addparameter(m, 'vaxS4', 439 - 1); % (March 14 = tv3 + 14)
    p(17) = addparameter(m, 'vaxJJ', 453 - 1); % enter J&J vaccine (March 26 = tvj + 28)
    p(18) = addparameter(m, 'vaxS5', 465 - 1); % (April 9 = tv3 + 32)
    p(19) = addparameter(m, 'vaxS6', 467 - 1); % (April 11 = tv4 + 14)
    p(20) = addparameter(m, 'vaxS7', 485 - 1); % (April 29 = tv4 + 32)
    p(21) = addparameter(m, 'vaxOff', 502 - 1); % end vaccinations (May 16 = tvacend)
    p(22) = addparameter(m, 'tsh', 213 - 1); % switch date for hospitalization l.o.s. parameters (July 31)
    p(23) = addparameter(m, 'tscc1', 170 - 1); % first switch date for cc parameters (June 18)
    p(24) = addparameter(m, 'tscc2', 273 - 1); % second switch date for cc parameters (September 29)
    p(25) = addparameter(m, 'vDuration', 1/730); % 1/duration of immunity from vaccine (1/vd)
    p(26) = addparameter(m, 'dimmuneA', 1/190); % 1/duration of immunity for asymptomatic recovery
    p(27) = addparameter(m, 'dimmuneI', 1/365); % 1/duration of immunity for symptomatic recovery
    
    % Time Dependent Parameters
    i = 30;
    bInd = zeros(1, Nr);
    NtInd = zeros(1, Nr);
    etaInd = zeros(Nr, Nr);
    thetaPiece = zeros(1, Nr);
    thetaInd = zeros(1, Nr);
    for j = regVec
        p(i) = addparameter(m, ['betaR' num2str(j)], 0);
        set(p(i), 'ConstantValue', false);
        bInd(j) = i;
        p(i + 1) = addparameter(m, ['bOptDay' num2str(tVec(1) + 1) 'inR' num2str(j)], 0.01);
        i = i + 2;
        p(i) = addparameter(m, ['Ntilde' num2str(j)], popByAge(j, 5));
        set(p(i), 'ConstantValue', false);
        NtInd(j) = i;
        i = i + 1;
        p(i) = addparameter(m, ['N' num2str(j)], popByAge(j, 5));
        i = i + 1;
        for k = regVec
            if j ~= k
                p(i) = addparameter(m, ['eta' num2str(j) 'to' num2str(k)], M0*contactMat(j, k, 1));
                set(p(i), 'ConstantValue', false);
                etaInd(j, k) = i;
                thetaPiece(j) = thetaPiece(j) + M0*contactMat(j, k, 1);
                i = i + 1;
            end
        end
        p(i) = addparameter(m, ['theta' num2str(j)], 1/(1 + thetaPiece(j)));
        set(p(i), 'ConstantValue', false);
        thetaInd(j) = i;
        i = i + 1;
    end
    
    for j = regVec % beta for j and j in j
        p(i) = addparameter(m, ['betaI' num2str(j) 'and' num2str(j) 'in' num2str(j)], p(bInd(j)).Value*p(2).Value*(p(thetaInd(j)).Value^2)/p(NtInd(j)).Value);
        set(p(i), 'ConstantValue', false);
        i = i + 1;
        p(i) = addparameter(m, ['betaA' num2str(j) 'and' num2str(j) 'in' num2str(j)], p(bInd(j)).Value*(p(thetaInd(j)).Value^2)/p(NtInd(j)).Value);
        set(p(i), 'ConstantValue', false);
        i = i + 1;
        g = regVec; g(g == j) = [];
        for k = g % beta for j and k in j
            p(i) = addparameter(m, ['betaI' num2str(j) 'and' num2str(k) 'in' num2str(j)], p(bInd(j)).Value*p(2).Value*p(thetaInd(j)).Value*p(thetaInd(k)).Value*p(etaInd(k, j)).Value/p(NtInd(j)).Value);
            set(p(i), 'ConstantValue', false);
            i = i + 1;
            p(i) = addparameter(m, ['betaA' num2str(j) 'and' num2str(k) 'in' num2str(j)], p(bInd(j)).Value*p(thetaInd(j)).Value*p(thetaInd(k)).Value*p(etaInd(k, j)).Value/p(NtInd(j)).Value);
            set(p(i), 'ConstantValue', false);
            i = i + 1;
            for jj = regVec 
                if jj == k
                    p(i) = addparameter(m, ['betaI' num2str(j) 'and' num2str(jj) 'in' num2str(k)], p(bInd(k)).Value*p(2).Value*p(thetaInd(j)).Value*p(thetaInd(jj)).Value*p(etaInd(j, k)).Value/p(NtInd(k)).Value);
                    set(p(i), 'ConstantValue', false);
                    i = i + 1;
                    p(i) = addparameter(m, ['betaA' num2str(j) 'and' num2str(jj) 'in' num2str(k)], p(bInd(k)).Value*p(thetaInd(j)).Value*p(thetaInd(jj)).Value*p(etaInd(j, k)).Value/p(NtInd(k)).Value);
                    set(p(i), 'ConstantValue', false);
                    i = i + 1;
                else
                    p(i) = addparameter(m, ['betaI' num2str(j) 'and' num2str(jj) 'in' num2str(k)], p(bInd(k)).Value*p(2).Value*p(thetaInd(j)).Value*p(thetaInd(jj)).Value*p(etaInd(j, k)).Value*p(etaInd(jj, k)).Value/p(NtInd(k)).Value);
                    set(p(i), 'ConstantValue', false);
                    i = i + 1;
                    p(i) = addparameter(m, ['betaA' num2str(j) 'and' num2str(jj) 'in' num2str(k)], p(bInd(k)).Value*p(thetaInd(j)).Value*p(thetaInd(jj)).Value*p(etaInd(j, k)).Value*p(etaInd(jj, k)).Value/p(NtInd(k)).Value);
                    set(p(i), 'ConstantValue', false);
                    i = i + 1;
                end
            end
        end
    end
    
    % Age specific parameters
    pS = [0.110023 0.35705 0.561205 0.774879]; % fraction of symptomatic cases
    dnh = [0.000013 0.0000822 0.000622705 0.027971448]; % death fraction for non-hospitalized patients % [0 0 0.00055021 0.021394]; 
    chosp = [0.0227054 0.03656527 0.050782 0.0789884]; % combined fraction hospitalized (hosp + cc)
    chlos = [5.8303 6.366447985 10.54149451 10.53665902]; % combined hospitalization length of stay %? (hosp.*hlos + cc.*clos)./chosp;
    dch = [0.005504587 0.014716188 0.057513214 0.160570555]; % combined death fraction for hospitalized patients %? (hosp.*dh + cc.*dc)./chosp;

    % Hospitalization parameter values after switch
    chlosb = [5.5747 5.231230631 8.486187615 8.136581647]; % B: July 31, tsh
    chospb = [0.0227054 0.03747127 0.049972 0.0724784]; % B: June 18, tscc1
    chospc = [0.0227054 0.02341227 0.03312364 0.10043694]; % C: September 29, tscc2
    
    j = 3000;
    for i = 1:Na % p(3) = alpha, p(4) = gamma
        p(j) = addparameter(m, ['iRecovery' num2str(i)], (p(4).Value)*(1 - chosp(i) - dnh(i)));
        set(p(j), 'ConstantValue', false);
        p(j + 1) = addparameter(m, ['iRecoveryB' num2str(i)], (p(4).Value)*(1 - chospb(i) - dnh(i)));
        p(j + 2) = addparameter(m, ['iRecoveryC' num2str(i)], (p(4).Value)*(1 - chospc(i) - dnh(i)));
        
        p(j + 3) = addparameter(m, ['hRecov' num2str(i)], (1 - dch(i))/chlos(i));
        set(p(j + 3), 'ConstantValue', false);
        p(j + 4) = addparameter(m, ['hRecovB' num2str(i)], (1 - dch(i))/chlosb(i));
        
        p(j + 5) = addparameter(m, ['iDeath' num2str(i)], dnh(i)*p(4).Value);
        p(j + 6) = addparameter(m, ['hDeath' num2str(i)], dch(i)/chlos(i));
        set(p(j + 6), 'ConstantValue', false);
        p(j + 7) = addparameter(m, ['hDeathB' num2str(i)], dch(i)/chlosb(i));
        
        p(j + 8) = addparameter(m, ['hHospitalized' num2str(i)], chosp(i)*p(4).Value);
        set(p(j + 8), 'ConstantValue', false);
        p(j + 9) = addparameter(m, ['hHospitalizedB' num2str(i)], chospb(i)*p(4).Value);
        p(j + 10) = addparameter(m, ['hHospitalizedC' num2str(i)], chospc(i)*p(4).Value);
        
        p(j + 11) = addparameter(m, ['eInfectious' num2str(i)], pS(i)/p(3).Value);
        p(j + 12) = addparameter(m, ['eAsymptomatic' num2str(i)], (1 - pS(i))/p(3).Value);
        j = j + 13;
    end
    
    % Vaccinations
    fc = [0.52 0.38 0.72]; % percent efficacy
    j = 3100;
    for i = 1:Na
        for k = regVec
            p(j) = addparameter(m, ['mrna1D' num2str(k) 'A' num2str(i)], 0);
            set(p(j), 'ConstantValue', false);
            p(j + 1) = addparameter(m, ['mrna2D' num2str(k) 'A' num2str(i)], 0);
            set(p(j + 1), 'ConstantValue', false);
            p(j + 2) = addparameter(m, ['jandj' num2str(k) 'A' num2str(i)], 0);
            set(p(j + 2), 'ConstantValue', false);
            p(j + 3) = addparameter(m, ['NR' num2str(k) 'A' num2str(i)], popByAge(k, i));
            p(j + 4) = addparameter(m, ['Nbar' num2str(k) 'A' num2str(i)], popByAge(k, i));
            set(p(j + 4), 'ConstantValue', false);
            p(j + 5) = addparameter(m, ['vaccRate' num2str(k) 'A' num2str(i)], (p(j).Value + p(j + 1).Value + p(j + 2).Value)/p(j + 4).Value);
            set(p(j + 5), 'ConstantValue', false);
            p(j + 6) = addparameter(m, ['mrna1DJan' num2str(k) 'A' num2str(i)], fc(1)*vacMat(1, i, k));
            p(j + 7) = addparameter(m, ['mrna1DFeb' num2str(k) 'A' num2str(i)], fc(1)*vacMat(2, i, k));
            p(j + 8) = addparameter(m, ['mrna1DMar' num2str(k) 'A' num2str(i)], fc(1)*vacMat(3, i, k));
            p(j + 9) = addparameter(m, ['mrna1DApr' num2str(k) 'A' num2str(i)], fc(1)*vacMat(4, i, k)); 
            p(j + 10) = addparameter(m, ['mrna2DJan' num2str(k) 'A' num2str(i)], fc(2)*vacMat(1, i, k));
            p(j + 11) = addparameter(m, ['mrna2DFeb' num2str(k) 'A' num2str(i)], fc(2)*vacMat(2, i, k));
            p(j + 12) = addparameter(m, ['mrna2DMar' num2str(k) 'A' num2str(i)], fc(2)*vacMat(3, i, k));
            p(j + 13) = addparameter(m, ['mrna2DApr' num2str(k) 'A' num2str(i)], fc(2)*vacMat(4, i, k));
            p(j + 14) = addparameter(m, ['jandjMar' num2str(k) 'A' num2str(i)], fc(3)*vacMat(5, i, k));
%             p(j + 15) = addparameter(m, ['jandjApr' num2str(k) 'A' num2str(i)], fc(3)*vacMat(6, i, k));
            j = j + 16;
        end
    end
    
    % Parameters for Optimization
    j = 3500;
    for i = 1:length(tVec)
        thetaPiece = zeros(1, Nr);
        thetaInd = zeros(1, Nr);
        etaPiece = zeros(Nr, Nr);
        p(j) = addparameter(m, ['tOptDay' num2str(tVec(i) + 1)], tVec(i));
        p(j + 1) = addparameter(m, ['tOpt' num2str(tVec(i) + 1) 'bd'], tVec(i) + evtBd);
        j = j + 2;
        for k = regVec
            
            for l = regVec
                if (l ~= k)
                    p(j) = addparameter(m, ['eta' num2str(k) 'to' num2str(l) 'Day' num2str(tVec(i) + 1)], M0*contactMat(k, l, i+1));
                    thetaPiece(k) = thetaPiece(k) + p(j).Value;
                    etaPiece(k, l) = p(NtInd(k)).Value*p(j).Value;
                    j = j + 1;
                end
            end
            
            p(j) = addparameter(m, ['theta' num2str(k) 'Day' num2str(tVec(i) + 1)], 1/(1 + thetaPiece(k)));
            etaPiece(k, :) = etaPiece(k, :)*p(j).Value;
            thetaInd(k) = j;
            j = j + 1;
            
        end
        
        for k = regVec
            p(j) = addparameter(m, ['Ntilde' num2str(k) 'Day' num2str(tVec(i) + 1)], p(thetaInd(k)).Value*p(NtInd(k)).Value + sum(etaPiece(:, k)));
            j = j + 1;
        end
        
    end
    
    for i = 2:length(tVec)
        for k = regVec
            p(j) = addparameter(m, ['bOptDay' num2str(tVec(i) + 1) 'inR' num2str(k)], bMatrix(k, (i - 1)));
            j = j + 1;
        end
    end

    
    %% Species and Compartments
    
    for i = regVec
        for j = 1:Na
            CO{i, j} = addcompartment(m, [ageCat{j} 'fromR' num2str(i)]);
            S{i, j} = addspecies(CO{i, j}, 'S');   % Susceptible
            E{i, j} = addspecies(CO{i, j}, 'E');   % Exposed
            I{i, j} = addspecies(CO{i, j}, 'I');   % (Symptomatic) Infectious  
            A{i, j} = addspecies(CO{i, j}, 'A');   % Asymptomatic Infectious
            RI{i, j} = addspecies(CO{i, j}, 'RI'); % Recovered (Symptomatic)
            RA{i, j} = addspecies(CO{i, j}, 'RA'); % Recovered (Asymptomatic)
            V{i, j} = addspecies(CO{i, j}, 'V');   % Vaccinated
            H{i, j} = addspecies(CO{i, j}, 'H');   % Hospitalized (all)
            D{i, j} = addspecies(CO{i, j}, 'D');   % Deceased
        end
        hosptot = addspecies(CO{i, 1}, 'hosptot'); 
    end
    htotTOT = addspecies(CO{1, 1}, 'htotTOT');
    htot65plus = addspecies(CO{1, 1}, 'htot65plus'); % (sometimes) used as an additional constraint when fitting lpha model
    
    %% Reactions
    
    for i = regVec
        for j = 1:Na
             
            % H Recovery
            H2RI{i, j} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.H -> ' ageCat{j} 'fromR' num2str(i) '.RI']);
            H2RI{i, j}.Name = [ageCat{j} 'fromR' num2str(i) ':H2RI'];
            H2RI{i, j}.addkineticlaw('MassAction');
            setparameter(H2RI{i, j}.KineticLaw, 'Forward Rate Parameter', ['hRecov' num2str(j)]);
            
            % Hospitalization Death
            H2D{i, j} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.H -> ' ageCat{j} 'fromR' num2str(i) '.D']);
            H2D{i, j}.Name = [ageCat{j} 'fromR' num2str(i) ':H2D'];
            H2D{i, j}.addkineticlaw('MassAction');
            setparameter(H2D{i, j}.KineticLaw, 'Forward Rate Parameter', ['hDeath' num2str(j)]);
            
            % Susceptible to Vaccinated
            S2V{i, j} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.S' ' -> ' ageCat{j} 'fromR' num2str(i) '.V']);
            S2V{i, j}.Name = [ageCat{j} 'fromR' num2str(i) ': S2V'];
            S2V{i, j}.addkineticlaw('MassAction');
            setparameter(S2V{i, j}.KineticLaw, 'Forward Rate Parameter', ['vaccRate' num2str(i) 'A' num2str(j)]);
            
            % Symptomatic Recovered to Vaccinated
            RI2V{i, j} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.RI' ' -> ' ageCat{j} 'fromR' num2str(i) '.V']);
            RI2V{i, j}.Name = [ageCat{j} 'fromR' num2str(i) ': RI2V'];
            RI2V{i, j}.addkineticlaw('MassAction');
            setparameter(RI2V{i, j}.KineticLaw, 'Forward Rate Parameter', ['vaccRate' num2str(i) 'A' num2str(j)]);
            
            % Asymptomatic Recovered to Vaccinated
            RA2V{i, j} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.RA' ' -> ' ageCat{j} 'fromR' num2str(i) '.V']);
            RA2V{i, j}.Name = [ageCat{j} 'fromR' num2str(i) ': RA2V'];
            RA2V{i, j}.addkineticlaw('MassAction');
            setparameter(RA2V{i, j}.KineticLaw, 'Forward Rate Parameter', ['vaccRate' num2str(i) 'A' num2str(j)]);

            % Symptomatic Recovered to Susceptible
            RI2S{i, j} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.RI -> ' ageCat{j} 'fromR' num2str(i) '.S']);
            RI2S{i, j}.Name = [ageCat{j} 'fromR' num2str(i) ': RI2S'];
            RI2S{i, j}.addkineticlaw('MassAction');
            setparameter(RI2S{i, j}.KineticLaw, 'Forward Rate Parameter', 'dimmuneI');

            % Asymptomatic Recovered to Susceptible
            RA2S{i, j} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.RA -> ' ageCat{j} 'fromR' num2str(i) '.S']);
            RA2S{i, j}.Name = [ageCat{j} 'fromR' num2str(i) ': RA2S'];
            RA2S{i, j}.addkineticlaw('MassAction');
            setparameter(RA2S{i, j}.KineticLaw, 'Forward Rate Parameter', 'dimmuneA');

            % Vaccinated to Susceptible
            V2S{i, j} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.V -> ' ageCat{j} 'fromR' num2str(i) '.S']);
            V2S{i, j}.Name = [ageCat{j} 'fromR' num2str(i) ': V2S'];
            V2S{i, j}.addkineticlaw('MassAction');
            setparameter(V2S{i, j}.KineticLaw, 'Forward Rate Parameter', 'vDuration');
            
            % Asymptomatic Recovery
            A2R{i, j} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.A -> ' ageCat{j} 'fromR' num2str(i) '.RA']);
            A2R{i, j}.Name = [ageCat{j} 'fromR' num2str(i) ': A2R'];
            A2R{i, j}.addkineticlaw('MassAction');
            setparameter(A2R{i, j}.KineticLaw, 'Forward Rate Parameter', 'gamma');

            % (Symptomatic) Infectious Recovery
            I2R{i, j} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.I -> ' ageCat{j} 'fromR' num2str(i) '.RI']);
            I2R{i, j}.Name = [ageCat{j} 'fromR' num2str(i) ': I2R'];
            I2R{i, j}.addkineticlaw('MassAction');
            setparameter(I2R{i, j}.KineticLaw, 'Forward Rate Parameter', ['iRecovery' num2str(j)]);

            % Non-Hospitalized Death
            I2D{i, j} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.I -> ' ageCat{j} 'fromR' num2str(i) '.D']);
            I2D{i, j}.Name = [ageCat{j} 'fromR' num2str(i) ': I2D'];
            I2D{i, j}.addkineticlaw('MassAction');
            setparameter(I2D{i, j}.KineticLaw, 'Forward Rate Parameter', ['iDeath' num2str(j)]);

            % Exposed to (Symptomatic) Infectious
            E2I{i, j} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.E -> ' ageCat{j} 'fromR' num2str(i) '.I']);
            E2I{i, j}.Name = [ageCat{j} 'fromR' num2str(i) ': E2I'];
            E2I{i, j}.addkineticlaw('MassAction');
            setparameter(E2I{i, j}.KineticLaw, 'Forward Rate Parameter', ['eInfectious' num2str(j)]);

            % Exposed to Asymptomatic Infectious
            E2A{i, j} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.E -> ' ageCat{j} 'fromR' num2str(i) '.A']);
            E2A{i, j}.Name = [ageCat{j} 'fromR' num2str(i) ': E2A'];
            E2A{i, j}.addkineticlaw('MassAction');
            setparameter(E2A{i, j}.KineticLaw, 'Forward Rate Parameter', ['eAsymptomatic' num2str(j)]);

            % I2H
            I2H{i, j} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.I -> ' ageCat{j} 'fromR' num2str(i) '.H']);
            I2H{i, j}.Name = [ageCat{j} 'fromR' num2str(i) ':I2H'];
            I2H{i, j}.addkineticlaw('MassAction');
            setparameter(I2H{i, j}.KineticLaw, 'Forward Rate Parameter', ['hHospitalized' num2str(j)]);

            % Nonlinear terms: Susceptible to Exposed
            for s = 1:Na

                % Exposed via Symptomatic Infectious Individual
                EbyI{j, s, i, i, i} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.S + ' ageCat{s} 'fromR' num2str(i) '.I -> ' ...
                    ageCat{s} 'fromR' num2str(i) '.I + ' ageCat{j} 'fromR' num2str(i) '.E']);
                EbyI{j, s, i, i, i}.Name = [ageCat{j} 'fromR' num2str(i) '.S - ' ageCat{s} 'fromR' num2str(i) '.I: Infection in R' num2str(i)];
                EbyI{j, s, i, i, i}.addkineticlaw('MassAction');
                setparameter(EbyI{j, s, i, i, i}.KineticLaw, 'Forward Rate Parameter', ['betaI' num2str(i) 'and' num2str(i) 'in' num2str(i)]);

                % Exposed via Asymptomatic Infectious Individual
                EbyA{j, s, i, i, i} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.S + ' ageCat{s} 'fromR' num2str(i) '.A -> ' ...
                    ageCat{s} 'fromR' num2str(i) '.A + ' ageCat{j} 'fromR' num2str(i) '.E']);
                EbyA{j, s, i, i, i}.Name = [ageCat{j} 'fromR' num2str(i) '.S - ' ageCat{s} 'fromR' num2str(i) '.A: Infection in R' num2str(i)];
                EbyA{j, s, i, i, i}.addkineticlaw('MassAction');
                setparameter(EbyA{j, s, i, i, i}.KineticLaw, 'Forward Rate Parameter', ['betaA' num2str(i) 'and' num2str(i) 'in' num2str(i)]);
                
                g = regVec; g(g == i) = [];
                for k = g
                    
                    % Exposed via Symptomatic Infectious Visitor to i from k
                    EbyI{j, s, i, k, i} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.S + ' ageCat{s} 'fromR' num2str(k) '.I -> ' ...
                        ageCat{s} 'fromR' num2str(k) '.I + ' ageCat{j} 'fromR' num2str(i) '.E']);
                    EbyI{j, s, i, k, i}.Name = [ageCat{j} 'fromR' num2str(i) '.S - ' ageCat{s} 'fromR' num2str(k) '.I: Infection in R' num2str(i)];
                    EbyI{j, s, i, k, i}.addkineticlaw('MassAction');
                    setparameter(EbyI{j, s, i, k, i}.KineticLaw, 'Forward Rate Parameter', ['betaI' num2str(i) 'and' num2str(k) 'in' num2str(i)]);

                    % Exposed via Asymptomatic Infectious Visitor to i from k
                    EbyA{j, s, i, k, i} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.S + ' ageCat{s} 'fromR' num2str(k) '.A -> ' ...
                        ageCat{s} 'fromR' num2str(k) '.A + ' ageCat{j} 'fromR' num2str(i) '.E']);
                    EbyA{j, s, i, k, i}.Name = [ageCat{j} 'fromR' num2str(i) '.S - ' ageCat{s} 'fromR' num2str(k) '.A: Infection in R' num2str(i)];
                    EbyA{j, s, i, k, i}.addkineticlaw('MassAction');
                    setparameter(EbyA{j, s, i, k, i}.KineticLaw, 'Forward Rate Parameter', ['betaA' num2str(i) 'and' num2str(k) 'in' num2str(i)]);
                    
                    for jj = regVec
                        
                        % Exposed via Symptomatic Infectious Visitor to k from jj
                        EbyI{j, s, i, jj, k} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.S + ' ageCat{s} 'fromR' num2str(jj) '.I -> ' ...
                            ageCat{s} 'fromR' num2str(jj) '.I + ' ageCat{j} 'fromR' num2str(i) '.E']);
                        EbyI{j, s, i, jj, k}.Name = [ageCat{j} 'fromR' num2str(i) '.S - ' ageCat{s} 'fromR' num2str(jj) '.I: Infection in R' num2str(k)];
                        EbyI{j, s, i, jj, k}.addkineticlaw('MassAction');
                        setparameter(EbyI{j, s, i, jj, k}.KineticLaw, 'Forward Rate Parameter', ['betaI' num2str(i) 'and' num2str(jj) 'in' num2str(k)]);

                        % Exposed via Asymptomatic Infectious Visitor to k from jj
                        EbyA{j, s, i, jj, k} = addreaction(m, [ageCat{j} 'fromR' num2str(i) '.S + ' ageCat{s} 'fromR' num2str(jj) '.A -> ' ...
                            ageCat{s} 'fromR' num2str(jj) '.A + ' ageCat{j} 'fromR' num2str(i) '.E']);
                        EbyA{j, s, i, jj, k}.Name = [ageCat{j} 'fromR' num2str(i) '.S - ' ageCat{s} 'fromR' num2str(jj) '.A: Infection in R' num2str(k)];
                        EbyA{j, s, i, jj, k}.addkineticlaw('MassAction');
                        setparameter(EbyA{j, s, i, jj, k}.KineticLaw, 'Forward Rate Parameter', ['betaA' num2str(i) 'and' num2str(jj) 'in' num2str(k)]);
                        
                    end
                    
                end

            end

        end
    end

    %% Initial Conditions 
    
    k = 1;
    for i = 1:Na:(Na*length(regVec))
        for j = 1:Na
            m.Compartment(i + (j - 1)).Species(1).InitialAmount = popByAge(k, j);
        end
        k = k + 1;
    end

    %% Events
    
    % Initial infection January 24, 2020
    for i = regVec %rI0
        evt1(i) = addevent(m, '(time >= ioBegin) && (time < ioBd)', {});
        set(evt1(i), 'Name', ['Initial infection (day 24 = 01/24/2020) in region ' num2str(i)]);
        set(evt1(i), 'EventFcns', {['[' ageCat{1} 'fromR' num2str(i) '].S = ([' ageCat{1} 'fromR' num2str(i) '].S - 1)'], ...
            ['[' ageCat{1} 'fromR' num2str(i) '].I = I0']});
    end

    % Changes in hospitalization parameters
    for i = 1:Na
        evt2(i) = addevent(m, 'time >= tsh', {}); % July 31
        set(evt2(i), 'Name', 'Change in hospitalization length of stay');
        set(evt2(i), 'EventFcns', {['hRecov' num2str(i) ' = hRecovB' num2str(i)], ...
            ['hDeath' num2str(i) ' = hDeathB' num2str(i)]});
    end
    
    for i = 1:Na
        evt3(i) = addevent(m, 'time >= tscc1', {}); % June 18
        set(evt3(i), 'Name', 'First change in fraction hospitalized');
        set(evt3(i), 'EventFcns', {['iRecovery' num2str(i) ' = iRecoveryB' num2str(i)], ...
            ['hHospitalized' num2str(i) ' = hHospitalizedB' num2str(i)]});
    end
    
    for i = 1:Na
        evt4(i) = addevent(m, 'time >= tscc2', {}); % September 29
        set(evt4(i), 'Name', 'Second change in fraction hospitalized');
        set(evt4(i), 'EventFcns', {['iRecovery' num2str(i) ' = iRecoveryC' num2str(i)], ...
            ['hHospitalized' num2str(i) ' = hHospitalizedC' num2str(i)]});
    end
    
    % Vaccination beginning, and changes in distribution
%     vTimes = {'vaxOn', 'vaxS1', 'vaxS2', 'vaxS3', 'vaxS4', 'vaxJJ', 'vaxS5', 'vaxS6', 'vaxOff'};
    vTags = {'Jan', 'Feb', 'Mar', 'Apr'};
    vTimes = {'vaxOn', 'vaxS2', 'vaxS4', 'vaxS6'};
    eCounter = 1;
    for i = 1:numel(vTimes)
        for k = regVec
            evt5(eCounter) = addevent(m, ['time >= ' vTimes{i}], {});
            set(evt5(eCounter), 'Name', ['Change in mRna (first dose) vaccine distribution: region ' num2str(k) ', no. ' num2str(i)]);
            set(evt5(eCounter), 'EventFcns', {['mrna1D' num2str(k) 'A1 = mrna1D' vTags{i} num2str(k) 'A1'], ...
                ['mrna1D' num2str(k) 'A2 = mrna1D' vTags{i} num2str(k) 'A2'], ...
                ['mrna1D' num2str(k) 'A3 = mrna1D' vTags{i} num2str(k) 'A3'], ...
                ['mrna1D' num2str(k) 'A4 = mrna1D' vTags{i} num2str(k) 'A4']});
            eCounter = eCounter + 1;
        end
    end
    
    vTimes = {'vaxS1', 'vaxS3', 'vaxS5', 'vaxS7'};
    for i = 1:numel(vTimes)
        for k = regVec
            evt5(eCounter) = addevent(m, ['time >= ' vTimes{i}], {});
            set(evt5(eCounter), 'Name', ['Change in mRna (second dose) vaccine distribution: region ' num2str(k) ', no. ' num2str(i)]);
            set(evt5(eCounter), 'EventFcns', {['mrna2D' num2str(k) 'A1 = mrna2D' vTags{i} num2str(k) 'A1'], ...
                ['mrna2D' num2str(k) 'A2 = mrna2D' vTags{i} num2str(k) 'A2'], ...
                ['mrna2D' num2str(k) 'A3 = mrna2D' vTags{i} num2str(k) 'A3'], ...
                ['mrna2D' num2str(k) 'A4 = mrna2D' vTags{i} num2str(k) 'A4']});
            eCounter = eCounter + 1;
        end
    end
    
    vTimes = {'vaxJJ'};
    for i = 1:numel(vTimes)
        for k = regVec
            evt5(eCounter) = addevent(m, ['time >= ' vTimes{i}], {});
            set(evt5(eCounter), 'Name', ['Change in J&J vaccine distribution: region ' num2str(k) ', no. ' num2str(i)]);
            set(evt5(eCounter), 'EventFcns', {['jandj' num2str(k) 'A1 = jandj' vTags{i + 2} num2str(k) 'A1'], ...
                ['jandj' num2str(k) 'A2 = jandj' vTags{i + 2} num2str(k) 'A2'], ...
                ['jandj' num2str(k) 'A3 = jandj' vTags{i + 2} num2str(k) 'A3'], ...
                ['jandj' num2str(k) 'A4 = jandj' vTags{i + 2} num2str(k) 'A4']});
            eCounter = eCounter + 1;
        end
    end
    
    % Bi-weekly changes in mobility and betas
    eCounter = 1;
    for i = 1:length(tVec)
        eFunc = cell(1, 1);
        fCount = 1;
        for j = regVec
            eFunc{fCount} = ['betaR' num2str(j) ' = bOptDay' num2str(tVec(i) + 1) 'inR' num2str(j)];
            eFunc{fCount + 1} = ['theta' num2str(j) ' = theta' num2str(j) 'Day' num2str(tVec(i) + 1)];
            eFunc{fCount + 2} = ['Ntilde' num2str(j) ' = Ntilde' num2str(j) 'Day' num2str(tVec(i) + 1)];
            eFunc{fCount + 3} = ['betaI' num2str(j) 'and' num2str(j) 'in' num2str(j) ' = bOptDay' num2str(tVec(i) + 1) 'inR' num2str(j) '*lambda', ...
                '*theta' num2str(j) 'Day' num2str(tVec(i) + 1) '*theta' num2str(j) 'Day' num2str(tVec(i) + 1) '/Ntilde' num2str(j) 'Day' num2str(tVec(i) + 1)];
            eFunc{fCount + 4} = ['betaA' num2str(j) 'and' num2str(j) 'in' num2str(j) ' = bOptDay' num2str(tVec(i) + 1) 'inR' num2str(j), ...
                '*theta' num2str(j) 'Day' num2str(tVec(i) + 1) '*theta' num2str(j) 'Day' num2str(tVec(i) + 1) '/Ntilde' num2str(j) 'Day' num2str(tVec(i) + 1)];
            fCount = fCount + 5;
            g = regVec; g(g == j) = [];
            for k = g 
                eFunc{fCount} = ['eta' num2str(j) 'to' num2str(k) ' = eta' num2str(j) 'to' num2str(k) 'Day' num2str(tVec(i) + 1)];
                fCount = fCount + 1;
            end
        end
        
        for j = regVec
            g = regVec; g(g == j) = [];
            for k = g % beta for j and k in j
                eFunc{fCount} = ['betaI' num2str(j) 'and' num2str(k) 'in' num2str(j) ' = bOptDay' num2str(tVec(i) + 1) 'inR' num2str(j) '*lambda', ...
                    '*theta' num2str(j) 'Day' num2str(tVec(i) + 1) '*theta' num2str(k) 'Day' num2str(tVec(i) + 1) '*eta' num2str(k) 'to' num2str(j) 'Day' num2str(tVec(i) + 1), ...
                    '/Ntilde' num2str(j) 'Day' num2str(tVec(i) + 1)];
                fCount = fCount + 1;
                eFunc{fCount} = ['betaA' num2str(j) 'and' num2str(k) 'in' num2str(j) ' = bOptDay' num2str(tVec(i) + 1) 'inR' num2str(j), ...
                    '*theta' num2str(j) 'Day' num2str(tVec(i) + 1) '*theta' num2str(k) 'Day' num2str(tVec(i) + 1) '*eta' num2str(k) 'to' num2str(j) 'Day' num2str(tVec(i) + 1), ...
                    '/Ntilde' num2str(j) 'Day' num2str(tVec(i) + 1)];
                fCount = fCount + 1;
                for jj = regVec % beta for j and jj in k
                    if jj == k
                        eFunc{fCount} = ['betaI' num2str(j) 'and' num2str(jj) 'in' num2str(k) ' = bOptDay' num2str(tVec(i) + 1) 'inR' num2str(k) '*lambda', ...
                            '*theta' num2str(j) 'Day' num2str(tVec(i) + 1) '*theta' num2str(jj) 'Day' num2str(tVec(i) + 1) '*eta' num2str(j) 'to' num2str(k) 'Day' num2str(tVec(i) + 1), ...
                            '/Ntilde' num2str(j) 'Day' num2str(tVec(i) + 1)];
                        fCount = fCount + 1;
                        eFunc{fCount} = ['betaA' num2str(j) 'and' num2str(jj) 'in' num2str(k) ' = bOptDay' num2str(tVec(i) + 1) 'inR' num2str(k), ...
                            '*theta' num2str(j) 'Day' num2str(tVec(i) + 1) '*theta' num2str(jj) 'Day' num2str(tVec(i) + 1) '*eta' num2str(j) 'to' num2str(k) 'Day' num2str(tVec(i) + 1), ...
                            '/Ntilde' num2str(k) 'Day' num2str(tVec(i) + 1)];
                        fCount = fCount + 1;
                    else
                        eFunc{fCount} = ['betaI' num2str(j) 'and' num2str(jj) 'in' num2str(k) ' = bOptDay' num2str(tVec(i) + 1) 'inR' num2str(k) '*lambda', ...
                            '*theta' num2str(j) 'Day' num2str(tVec(i) + 1) '*theta' num2str(jj) 'Day' num2str(tVec(i) + 1) '*eta' num2str(j) 'to' num2str(k) 'Day' num2str(tVec(i) + 1), ...
                            '*eta' num2str(jj) 'to' num2str(k) 'Day' num2str(tVec(i) + 1) '/Ntilde' num2str(j) 'Day' num2str(tVec(i) + 1)];
                        fCount = fCount + 1;
                        eFunc{fCount} = ['betaA' num2str(j) 'and' num2str(jj) 'in' num2str(k) ' = bOptDay' num2str(tVec(i) + 1) 'inR' num2str(k), ...
                            '*theta' num2str(j) 'Day' num2str(tVec(i) + 1) '*theta' num2str(jj) 'Day' num2str(tVec(i) + 1) '*eta' num2str(j) 'to' num2str(k) 'Day' num2str(tVec(i) + 1), ...
                            '*eta' num2str(jj) 'to' num2str(k) 'Day' num2str(tVec(i) + 1) '/Ntilde' num2str(k) 'Day' num2str(tVec(i) + 1)];
                        fCount = fCount + 1;
                    end
                end
            end
        end

        evt7(eCounter) = addevent(m, ['(time >= tOptDay' num2str(tVec(i) + 1) ') && (time < tOpt' num2str(tVec(i) + 1) 'bd)'], {});
        set(evt7(eCounter), 'Name', ['Mobility and Beta Changes for Day ' num2str(tVec(i) + 1)]);
        set(evt7(eCounter), 'EventFcns', eFunc);
        eCounter = eCounter + 1;
        
    end 
    
    %% Fitting and Optimization
    
    strOpt1 = ''; strOpt5 = '';
    for k = regVec
        
        % Total Hospitalizations
        strOpt5 = [strOpt5 ' + [' ageCat{4} 'fromR' num2str(k) '].H'];
        strOpt2 = ['[' ageCat{1} 'fromR' num2str(k) '].H'];
        for i = 2:Na
            strOpt2 = [strOpt2 ' + [' ageCat{i} 'fromR' num2str(k) '].H'];
        end
        strOpt1 = [strOpt1 ' + ' strOpt2];
        strOpt2 = ['[' ageCat{1} 'fromR' num2str(k) '].hosptot = ' strOpt2];
        addrule(m, strOpt2, 'RuleType', 'repeatedAssignment');
        
        % Vaccination rates and Nbar (population eligible for vaccination)
        for i = 1:Na
            strOpt3 = ['[' ageCat{i} 'fromR' num2str(k) '].H + [' ageCat{i} 'fromR' num2str(k) '].D + [' ageCat{i} 'fromR' num2str(k) '].V'];
            strOpt3 = ['theta' num2str(k) '*NR' num2str(k) 'A' num2str(i) ' - (' strOpt3 ')'];
            strOpt3 = ['Nbar' num2str(k) 'A' num2str(i) ' = ' strOpt3];
            addrule(m, strOpt3, 'RuleType', 'repeatedAssignment');
            strOpt4 = ['vaccRate' num2str(k) 'A' num2str(i) ' = (mrna1D' num2str(k) 'A' num2str(i) ' + mrna2D' num2str(k) 'A' num2str(i) ' + jandj', ...
                num2str(k) 'A' num2str(i) ')/max(1, Nbar' num2str(k) 'A' num2str(i) ')'];
            addrule(m, strOpt4, 'RuleType', 'repeatedAssignment');
        end
        
    end

    strOpt1 = ['[' ageCat{1} 'fromR1].htotTOT = ' strOpt1];
    addrule(m, strOpt1, 'RuleType', 'repeatedAssignment'); 
    strOpt5 = ['[' ageCat{1} 'fromR1].htot65plus = ' strOpt5];
    addrule(m, strOpt5, 'RuleType', 'repeatedAssignment');      
    
end
