%%%%%%%%%%%%%%%%% DISTRIBUTIONS %%%%%%%%%%%%%%%%%%%%%%%%
x = 0:.1:15;
basic = x;
traddist = normpdf(x, 5, 1);

s1 = normpdf(x, 3, .5);
s2 = normpdf(x, 7, .5);
scombined = s1+s2;
socialdist = scombined / trapz(x,scombined);
inddist = normpdf(x, 5.5, 1);


%%%%%%%%%%%%%%%% END DISTRIBUTIONS %%%%%%%%%%%%%%%%%%%%%

misinfdist = normpdf(x, 9, 0.0000001);
BETA = 0.00000020205;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
months = 3;
days = months * 30;
t = (days)*10; %(days * 10) 
step_size = 0.1;  % Step size
total_length = round(t / step_size);

pts = 0.576;  % Some constant
pst = 2.302;  % Some constant
pti = 0.329;  % Some constant
psi = 0.329;  % Some constant

klostrad = .657;
tradperc = .5;
klostsocial = 1.53;
socialperc = .3;
klostind = .153;
indperc = .15;
klostmis = .461;
misperc = .05;



tradists = zeros(total_length, length(basic));
tradists(1, :) = traddist; 

socialdists = zeros(total_length, length(basic)); 
socialdists(1, :) = socialdist; 

inddists = zeros(total_length, length(basic));  
inddists(1, :) = inddist; 

totaldists = zeros(total_length, length(basic));  
totaldists(1, :) = traddist + socialdist + inddist;


x_axis = [0];  % Will store the days incremented by whole numbers

alltrad = [0.0000000001];  % Initialize traditional array
allsocial = [0.000000001];  % Initialize social array
allindividual = [0.0000000001];  % Initialize individual array
allmis = [0.00000000000001];


mistrad = [0.0000000000001];
missocial = [0.00000000001];
misindividual = [0.000000000001];

indret = [0.00000001];
tradret = [0];

current_x = 0;  % Current position in time
i = 1;  % Index for the arrays
total = [0];

bigX = 1.16;
kreti =.05003;
krettrad = kreti/bigX;
kretsoc = kreti/6;
kretmis = kreti/2;

tol = .005;  % Define a small tolerance
misgrowth = true;
miscounter = 0;

injectioncounter = -2;
%%%%%%%%%%%%%%%%%%%%%%% GENERATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% While loop for calculations
while current_x < t
    % Convert current time to days for x-axis
    if mod(i-1, (1/step_size)) == 0  % Every 100 steps, increment the day
        x_axis(end + 1) = (i-1) / (1/step_size);  % Increment day by 1
    end
    injecttime = 5000;
    injectamount = .05;



    % Temporary variables to store next values
    new_trad = alltrad(i);
    new_social = allsocial(i);
    new_individual = allindividual(i);
    new_mis = allmis(i);
    newmistrad = 0;
    newmissocial =0;
    newmisindividual = 0;
    if miscounter < 1
        if abs(mod(i*step_size-.2, 7)) < tol
            misgrowth = true;
        end
        if misgrowth == true
            new_mis = new_mis + (klostmis * misperc - klostmis * new_mis) * step_size;
            miscounter = miscounter + step_size;
        end
    else 
        miscounter = 0;
        misgrowth = false;
    end
    

    
    if abs(mod(i*step_size-.1, 8.1+(7*injectioncounter))) < tol*.1
    
        alltrad(i) = alltrad(i) + allmis(i)/2;
        allsocial(i) = allsocial(i) +  allmis(i)/2;
        newmistrad = mistrad(i) +  allmis(i)/2;
        newmissocial = missocial(i) + allmis(i)/2;
        allmis(i) = 0;
        new_mis = 0;
        injectioncounter = injectioncounter + 1;
    end

    newmistrad = newmistrad + mistrad(i);
    newmissocial = newmissocial + missocial(i);
    newmisindividual = newmisindividual + misindividual(i);

    newtraddist = tradists(i, :);
    newsocialdist = socialdists(i, :);
    newinddist = inddists(i, :);



    diffTS = alltrad(i) - allsocial(i);  % Difference between trad and social
    diffTI = alltrad(i) - allindividual(i);  % Difference between trad and individual
    diffIS = allindividual(i) - allsocial(i);  % Difference between individual and social

    if mod((i*step_size), 12*30) == 0 % @ 12 months
        scale = trapz(x, totaldists(i,:)) * BETA;

        newtraddist = tradists(i, :) + scale * misinfdist(1, :)/2 ;
        newsocialdist = socialdists(i,:) + scale * misinfdist(1, :)/2;

        new_trad = new_trad + BETA * total(i)/2;
        new_social = new_social + BETA * total(i)/2;

    end
    
    % Update for traditional vs social
    if diffTS > 0
        changeT =  (pts * diffTS) * step_size + (klostrad * tradperc - klostrad * new_trad) * step_size;
    
        changeFromFlow = (pts * diffTS) * step_size / alltrad(i);
        changeTfromGrowth = (klostrad * tradperc - klostrad * new_trad) * step_size / alltrad(i);
        changeS = (pts * diffTS) * step_size + (klostsocial * socialperc - klostsocial * new_social) * step_size;
        changeSfromGrowth = (klostsocial * socialperc - klostsocial * new_social) * step_size / allsocial(i);

        if  newmistrad - changeFromFlow * mistrad(i) > tol
            newmistrad = newmistrad - changeFromFlow * mistrad(i);
            newmissocial = newmissocial + changeFromFlow * mistrad(i);
        end
        

        newtraddist(1, :) = newtraddist(1, :) - changeFromFlow * tradists(i,:) + changeTfromGrowth * tradists(i,:);
        newsocialdist(1, :) = newsocialdist(1, :) + changeFromFlow * tradists(i, :) + changeSfromGrowth * socialdists(i, :);

        new_trad = new_trad - (pts * diffTS) * step_size + (klostrad * tradperc - klostrad * new_trad) * step_size;
        new_social = new_social + (pts * diffTS) * step_size + (klostsocial * socialperc - klostsocial * new_social) * step_size;

    elseif diffTS < 0
        changeT = (pst * -diffTS) * step_size + (klostrad * tradperc - klostrad * new_trad) * step_size;
        changeS =  (pst * -diffTS) * step_size + (klostsocial * socialperc -klostsocial * new_social) * step_size;

        changeFromFlow = (pst * -diffTS) / allsocial(i);
        changeTfromGrowth = (klostrad * tradperc - klostrad * new_trad) * step_size / alltrad(i);
        changeSfromGrowth = (klostsocial * socialperc -klostsocial * new_social) * step_size / allsocial(i);
        
        if newmissocial - changeFromFlow*allsocial(i) > tol
            newmistrad = newmistrad + changeFromFlow*allsocial(i);
            newmissocial = newmissocial - changeFromFlow*allsocial(i);
        end

        newtraddist(1, :) = newtraddist(1, :) + changeFromFlow*socialdists(i, :) + changeTfromGrowth * tradists(i, :);
        newsocialdist(1, :) = newsocialdist(1, :) - changeFromFlow*socialdists(i, :) + changeSfromGrowth * socialdists(i, :);

        new_trad = new_trad + (pst * -diffTS) * step_size + (klostrad * tradperc - klostrad * new_trad) * step_size;
        new_social = new_social - (pst * -diffTS) * step_size + (klostsocial * socialperc -klostsocial * new_social) * step_size;

    elseif diffTS == 0
        changeT = (klostrad * tradperc - klostrad * new_trad) * step_size / alltrad(i);
        changeS = (klostsocial * socialperc - klostsocial * new_social) * step_size / allsocial(i);

        newtraddist(1, :) = newtraddist(1, :) + (changeT)*tradists(i, :);
        newsocialdist(1, :) = newsocialdist(1, :) + (changeS)*socialdists(i, :);

        new_trad = new_trad + (klostrad * tradperc - klostrad * new_trad) * step_size;
        new_social = new_social + (klostsocial * socialperc - klostsocial * new_social) * step_size;
    end

    % Update for individual vs social

    if diffIS < 0
        changeI = (klostind * indperc - klostind * new_individual) * step_size / allindividual(i);
        changeS = (klostsocial * socialperc - klostsocial * new_social) * step_size / allsocial(i);
        changeFromFlow = (psi * -diffIS) / allsocial(i);
        
        if  newmissocial - changeFromFlow * allsocial(i) > tol
            newmisindividual = newmisindividual + changeFromFlow * allsocial(i);
            newmissocial = newmissocial - changeFromFlow * allsocial(i);
        end

        newinddist(1, :) = newinddist(1, :) + changeFromFlow*socialdists(i, :) + changeI*(inddists(i,:));
        newsocialdist(1, :) = newsocialdist(1, :) - changeFromFlow*socialdists(i, :) + changeS*socialdists(i, :);

        new_individual = new_individual + (psi * -diffIS) * step_size + (klostind * indperc - klostind * new_individual) * step_size;
        new_social = new_social - (psi * -diffIS) * step_size + (klostsocial * socialperc - klostsocial * new_social) * step_size;
    elseif diffIS == 0
        changeI = (klostind * indperc - klostind * new_individual) * step_size / allindividual(i);
        changeS = (klostsocial * socialperc - klostsocial * new_social) * step_size / allsocial(i);

        newinddist(1, :) = newinddist(1, :) + (changeI)*inddists(i, :);
        newsocialdist(1, :) = newsocialdist(1, :) + (changeS)*socialdists(i, :); 

        new_individual = (klostind * indperc - klostind * new_individual) * step_size;
        new_social = (klostsocial * socialperc - klostsocial * new_social) * step_size;
    end
    
    % Update for traditional vs individual
    if diffTI > 0
        changeT =  (klostrad * tradperc -klostrad * new_trad) * step_size / alltrad(i);
        changeI =  (klostsocial * socialperc - klostsocial* new_individual) * step_size / allsocial(i);
        changeFromFlow = (pti * diffTI) * step_size / alltrad(i);

        newinddist(1, :) = newinddist(1, :) + changeFromFlow*tradists(i,:) + changeI*inddists(i,:);
        newtraddist(1, :) = newtraddist(1, :) - changeFromFlow*tradists(i,:) + changeT*tradists(i,:);
        
        if newmistrad - changeFromFlow * alltrad(i) > tol
            newmisindividual = newmisindividual + changeFromFlow * alltrad(i);
            newmistrad = newmistrad - changeFromFlow * alltrad(i);
        end

        new_trad = new_trad - (pti * diffTI) * step_size + (klostrad * tradperc -klostrad * new_trad) * step_size;
        new_individual = new_individual + (pti * diffTI) * step_size + (klostsocial * socialperc - klostsocial* new_individual) * step_size;
    elseif diffTI == 0
        changeI = (klostind * indperc - klostind * new_individual) * step_size / allindividual(i);
        changeT = (klostrad * tradperc - klostrad * new_trad) * step_size / alltrad(i);

        newinddist(1, :) = newinddist(1, :) + (changeI)*inddists(i, :);
        newtraddist(1, :) = newtraddist(1, :) + (changeT)*tradists(i, :);

        new_individual = (klostind * indperc - klostind * new_individual) * step_size;
        new_trad = (klostrad * tradperc - klostrad * new_trad) * step_size;
    end
    
    newmisindividual = newmisindividual + (- klostmis * newmisindividual) * step_size;
    new_mis = new_mis + (- klostmis * new_mis) * step_size;
    newmissocial = newmissocial + (- klostmis * newmissocial) * step_size;
    newmistrad = newmistrad + (- klostmis * newmistrad) * step_size;

   
    % Store the updated values in the arrays
    alltrad(i + 1) = new_trad;
    allsocial(i + 1) = new_social;
    allindividual(i + 1) = new_individual;
    %total(i+1) = new_trad + new_social + new_individual;

    mistrad(i+1) = newmistrad;
    missocial(i+1) =  newmissocial;
    misindividual(i+1) = newmisindividual;
    allmis(i+1) = new_mis;


    nettrad = new_trad - newmistrad;
    netsocial = new_social - newmissocial;
    netindividual = new_individual - newmisindividual;
    netmis = newmisindividual + newmissocial + newmistrad;

    total(i+1) = new_mis +  new_individual + new_social + new_trad;

    indret(i+1) = indret(i)+ (nettrad*krettrad + netsocial*kretsoc + netindividual*kreti + new_mis*kretmis) *step_size;

    tradret(i+1) = tradret(i) + nettrad*krettrad * step_size;

    tradists(i+1, :) = newtraddist(1, :);
    socialdists(i+1, :) = newsocialdist(1, :);
    inddists(i+1, :) = newinddist(1, :);
    totaldists(i+1, :) = newtraddist(1, :) + newsocialdist(1, :) + newinddist(1, :);

    % Update current time
    current_x = current_x + step_size;  
    i = i + 1; 
end

x_axis = x_axis * step_size;




% Plotting alltrad, allsocial, and allindividual vs x_axis
figure;
hold on;
plot(x_axis, alltrad(1:length(x_axis)), 'r', 'DisplayName', 'Traditional');
plot(x_axis, allsocial(1:length(x_axis)), 'b', 'DisplayName', 'Social');
plot(x_axis, allindividual(1:length(x_axis)), 'g', 'DisplayName', 'Individual');
plot(x_axis, total(1:length(x_axis)), 'DisplayName', 'Total');
% plot(x_axis, mistrad(1:length(x_axis)), ':', 'LineWidth', 1.5, 'DisplayName', "trad mis");
% plot(x_axis, missocial(1:length(x_axis)), ':', 'LineWidth', 1.5,'DisplayName', "social mis");
% plot(x_axis, misindividual(1:length(x_axis)), ':', 'LineWidth', 1.5,'DisplayName', "individual mis");
plot(x_axis, tradret(1:length(x_axis)), "DisplayName", "Trad Retained")
plot(x_axis, indret(1:length(x_axis)), 'LineWidth' , 1.5, 'DisplayName', 'retained');
plot(x_axis, allmis(1:length(x_axis)), 'LineWidth' , .1, 'DisplayName', 'injected MIS');

hold off;

% Adding titles and labels
title('INF Over Time');
xlabel('Time (days)');
ylabel('Proportion of INF expected');
legend show;
grid on;

disp(indret(90/step_size))
disp(tradret(90/step_size))
disp(tradret(90/step_size)/indret(90/step_size))