%%%%%%%%%%%%%%%%% DISTRIBUTIONS %%%%%%%%%%%%%%%%%%%%%%%%
x = 0:.5:15;
basic = x;
traddist = normpdf(x, 5, 1);

s1 = normpdf(x, 3, .5);
s2 = normpdf(x, 7, .5);
scombined = s1+s2;
socialdist = scombined / trapz(x,scombined);

p1 = .25 * normpdf(x, 3, 5);
p2 = .3 * normpdf(x, 7, .5);
p3 = .45 * normpdf(x, 5.5, 1.5);
combined = p1+p2+p3;
inddist = combined / trapz(x, combined);

%%%%%%%%%%%%%%%% END DISTRIBUTIONS %%%%%%%%%%%%%%%%%%%%%
ogmisinfdist = normpdf(x, 1, 0.001);
misinfdist = ogmisinfdist ;

BETA = 34;
ogbetaflag = true
BETAFLAG = ogbetaflag;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
months = 12;
days = months * 30;
t = (days+1)*10; %(days * 10) 
step_size = 0.1;  % Step size
total_length = round(t / step_size);

pts = 0.576;  % Some constant
pst = 2.302;  % Some constant
pti = 0.329;  % Some constant
psi = 0.329;  % Some constant

klostrad = .657;
tradperc = .45;
klostsocial = 1.53;
socialperc = .45;
klostind = .153;
indperc = .10;


tradists = zeros(total_length, length(basic));
tradists(1, :) = oldtraddist; 

socialdists = zeros(total_length, length(basic)); 
socialdists(1, :) = oldsocialdist; 

inddists = zeros(total_length, length(basic));  
inddists(1, :) = oldinddist; 

totaldists = zeros(total_length, length(basic));  
totaldists(1, :) = oldtraddist + oldsocialdist + oldinddist;

x_axis = [0];  % Will store the days incremented by whole numbers
alltrad = [oldtrad];  % Initialize traditional array
allsocial = [oldsocial];  % Initialize social array
allindividual = [oldindividual];  % Initialize individual array
indret = [oldret];
current_x = 0;  % Current position in time
i = 1;  % Index for the arrays
total = [oldtotal];

kreti = .9;
krettrad = .05;
kretsoc = .05;

kreti = .05003;
krettrad = kreti/ .9 * .05;
kretsoc = kreti/ .9 * .05;
%%%%%%%%%%%%%%%%%%%%%%% GENERATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while current_x < t
    % Convert current time to days for x-axis
    if mod(i-1, (1/step_size)) == 0  % Every 100 steps, increment the day
        x_axis(end + 1) = (i-1) / (1/step_size);  % Increment day by 1
    end

    diffTS = alltrad(i) - allsocial(i);  % Difference between trad and social
    diffTI = alltrad(i) - allindividual(i);  % Difference between trad and individual
    diffIS = allindividual(i) - allsocial(i);  % Difference between individual and social

    % Temporary variables to store next values
    new_trad = alltrad(i);
    new_social = allsocial(i);
    new_individual = allindividual(i);

    newtraddist = tradists(i, :);
    newsocialdist = socialdists(i, :);
    newinddist = inddists(i, :);


    new_social = new_social + (kgsocial - klostsocial * new_social)* socialperc *step_size ;
    new_trad = new_trad + (kgtrad  - klostrad * new_trad) * tradperc*step_size ;
    new_individual = new_individual + (kgind  - klostind * new_individual)* indperc*step_size ;

    changeI = (kgind  - klostind * new_individual)* indperc  *step_size/ allindividual(i);
    changeT = (kgtrad  - klostrad * new_trad)* tradperc *step_size / alltrad(i);
    changeS = (kgsocial  - klostsocial * new_social)* socialperc  * step_size/ allsocial(i);

    newtraddist = newtraddist + (changeT)*traddist;
    newinddist = newinddist + (changeI)*newinddist;
    newsocialdist = newsocialdist + (changeS)*socialdist;

    %Misinfinjection
    if mod((i*step_size), 12*30) == 0 && BETAFLAG == true % @ 12 months
        misinfdist= misinfdist * BETA * .01 * total(i);
        newtraddist = newtraddist + misinfdist/2 ;
        newsocialdist = newsocialdist  + misinfdist/2;
        new_trad = new_trad + BETA  * .01* total(i)/2;
        new_social = new_social + BETA  * .01* total(i)/2;
        
        BETAFLAG = ogbetaflag;
    end
 

    % Update for traditional vs social
    if diffTS > 0
        changeFromFlow = (pts * diffTS)  * step_size / alltrad(i);

            newtraddist = newtraddist - changeFromFlow * tradists(i,:);
            newsocialdist = newsocialdist + changeFromFlow * tradists(i, :);
    
            new_trad = new_trad - (pts * diffTS) * step_size;
            new_social = new_social + (pts * diffTS) * step_size;

    elseif diffTS < 0
        changeFromFlow = (pst * -diffTS) * step_size  / allsocial(i);

            newtraddist = newtraddist + changeFromFlow*socialdists(i, :);
            newsocialdist = newsocialdist - changeFromFlow*socialdists(i, :);
    
            new_trad = new_trad + (pst * -diffTS) * step_size;
            new_social = new_social - (pst * -diffTS) * step_size ;
    end
    
    if diffIS < 0
        changeFromFlow = (psi * -diffIS)* step_size  / allsocial(i);

        newinddist = newinddist + changeFromFlow*socialdists(i, :);
        newsocialdist = newsocialdist - changeFromFlow*socialdists(i, :);

        new_individual = new_individual + (psi * -diffIS) * step_size;
        new_social = new_social - (psi * -diffIS) * step_size;
    end

    if diffTI > 0
        changeFromFlow = (pti * diffTI) * step_size  / alltrad(i);

        newinddist = newinddist + changeFromFlow*tradists(i,:) ;
        newtraddist = newtraddist - changeFromFlow*tradists(i,:) ;
       
        new_trad = new_trad - (pti * diffTI)* step_size;
        new_individual = new_individual + (pti * diffTI)* step_size ;
    end

    alltrad(i + 1) = new_trad;
    allsocial(i + 1) = new_social;
    allindividual(i + 1) = new_individual;
    total(i+1) = new_trad + new_social + new_individual;
    
    indret(i+1) = indret(i)+ (new_trad * krettrad + new_social * kretsoc + new_individual * kreti)*step_size;
    
    tradists(i+1, :) = newtraddist(1, :);
    socialdists(i+1, :) = newsocialdist(1, :);
    inddists(i+1, :) = newinddist(1, :);

    totaldists(i+1, :) = newtraddist(1, :) + newsocialdist(1, :) + newinddist(1, :);

    current_x = current_x + step_size; 
    i = i + 1;  
end
x_axis = x_axis * step_size;

Atrad = alltrad;
Asocial = allsocial;
Aind = allindividual;
Atot = total;
Aret = indret;
POPA = inddists;


%%%%%%%%%%%%%%%%% POP B DISTRIBUTIONS %%%%%%%%%%%%%%%%%%%%%%%%
x = 0:.5:15;
basic = x;
traddist = normpdf(x, 5, 1);

socialdist = scombined / trapz(x,scombined);

p1 = .55 * normpdf(x, 4, 1);
p2 = .45 * normpdf(x, 6, 1);
combined = p1+p2;
inddist = combined / trapz(x, combined);

%%%%%%%%%%%%%%%% END DISTRIBUTIONS %%%%%%%%%%%%%%%%%%%%%

misinfdist = ogmisinfdist ;
BETAFLAG =true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
months = 12;
days = months * 30;
t = (days+1)*10; %(days * 10) 
step_size = 0.1;  % Step size
total_length = round(t / step_size);

pts = 0.576;  % Some constant
pst = 2.302;  % Some constant
pti = 0.329;  % Some constant
psi = 0.329;  % Some constant

klostrad = .657;
tradperc = .45;
klostsocial = 1.53;
socialperc = .45;
klostind = .153;
indperc = .10;

tradists = zeros(total_length, length(basic));
tradists(1, :) = oldtraddist; 

socialdists = zeros(total_length, length(basic)); 
socialdists(1, :) = oldsocialdist; 

inddists = zeros(total_length, length(basic));  
inddists(1, :) = oldinddist; 

totaldists = zeros(total_length, length(basic));  
totaldists(1, :) = oldtraddist + oldsocialdist + oldinddist;

x_axis = [0];  % Will store the days incremented by whole numbers
alltrad = [oldtrad];  % Initialize traditional array
allsocial = [oldsocial];  % Initialize social array
allindividual = [oldindividual];  % Initialize individual array
indret = [oldret];
current_x = 0;  % Current position in time
i = 1;  % Index for the arrays
total = [oldtotal];

kreti = .6;
krettrad = .25;
kretsoc = .15;


kreti = .05003;
krettrad = kreti/ .6 * .25;
kretsoc = kreti/ .6 * .15;
%%%%%%%%%%%%%%%%%%%%%%% GENERATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% While loop for calculations
while current_x < t
    % Convert current time to days for x-axis
    if mod(i-1, (1/step_size)) == 0  % Every 100 steps, increment the day
        x_axis(end + 1) = (i-1) / (1/step_size);  % Increment day by 1
    end

    diffTS = alltrad(i) - allsocial(i);  % Difference between trad and social
    diffTI = alltrad(i) - allindividual(i);  % Difference between trad and individual
    diffIS = allindividual(i) - allsocial(i);  % Difference between individual and social

    % Temporary variables to store next values
    new_trad = alltrad(i);
    new_social = allsocial(i);
    new_individual = allindividual(i);

    newtraddist = tradists(i, :);
    newsocialdist = socialdists(i, :);
    newinddist = inddists(i, :);


    new_social = new_social + (kgsocial - klostsocial * new_social)* socialperc *step_size ;
    new_trad = new_trad + (kgtrad  - klostrad * new_trad) * tradperc*step_size ;
    new_individual = new_individual + (kgind  - klostind * new_individual)* indperc*step_size ;

    changeI = (kgind  - klostind * new_individual)* indperc  *step_size/ allindividual(i);
    changeT = (kgtrad  - klostrad * new_trad)* tradperc *step_size / alltrad(i);
    changeS = (kgsocial  - klostsocial * new_social)* socialperc  * step_size/ allsocial(i);

    newtraddist = newtraddist + (changeT)*traddist;
    newinddist = newinddist + (changeI)*newinddist;
    newsocialdist = newsocialdist + (changeS)*socialdist;

    %Misinfinjection
    if mod((i*step_size), 12*30) == 0 && BETAFLAG == true % @ 12 months
        misinfdist= misinfdist * BETA * .01 * total(i);
        newtraddist = newtraddist + misinfdist/2 ;
        newsocialdist = newsocialdist  + misinfdist/2;
        new_trad = new_trad + BETA  * .01* total(i)/2;
        new_social = new_social + BETA  * .01* total(i)/2;
        
        BETAFLAG = false;
    end
 

    % Update for traditional vs social
    if diffTS > 0
        changeFromFlow = (pts * diffTS)  * step_size / alltrad(i);

            newtraddist = newtraddist - changeFromFlow * tradists(i,:);
            newsocialdist = newsocialdist + changeFromFlow * tradists(i, :);
    
            new_trad = new_trad - (pts * diffTS) * step_size;
            new_social = new_social + (pts * diffTS) * step_size;

    elseif diffTS < 0
        changeFromFlow = (pst * -diffTS) * step_size  / allsocial(i);

            newtraddist = newtraddist + changeFromFlow*socialdists(i, :);
            newsocialdist = newsocialdist - changeFromFlow*socialdists(i, :);
    
            new_trad = new_trad + (pst * -diffTS) * step_size;
            new_social = new_social - (pst * -diffTS) * step_size ;
    end
    
    if diffIS < 0
        changeFromFlow = (psi * -diffIS)* step_size  / allsocial(i);

        newinddist = newinddist + changeFromFlow*socialdists(i, :);
        newsocialdist = newsocialdist - changeFromFlow*socialdists(i, :);

        new_individual = new_individual + (psi * -diffIS) * step_size;
        new_social = new_social - (psi * -diffIS) * step_size;
    end

    if diffTI > 0
        changeFromFlow = (pti * diffTI) * step_size  / alltrad(i);

        newinddist = newinddist + changeFromFlow*tradists(i,:) ;
        newtraddist = newtraddist - changeFromFlow*tradists(i,:) ;
       
        new_trad = new_trad - (pti * diffTI)* step_size;
        new_individual = new_individual + (pti * diffTI)* step_size ;
    end

    alltrad(i + 1) = new_trad;
    allsocial(i + 1) = new_social;
    allindividual(i + 1) = new_individual;
    total(i+1) = new_trad + new_social + new_individual;
    
    indret(i+1) = indret(i)+ (new_trad * krettrad + new_social * kretsoc + new_individual * kreti)*step_size;
    
    tradists(i+1, :) = newtraddist(1, :);
    socialdists(i+1, :) = newsocialdist(1, :);
    inddists(i+1, :) = newinddist(1, :);

    totaldists(i+1, :) = newtraddist(1, :) + newsocialdist(1, :) + newinddist(1, :);

    current_x = current_x + step_size; 
    i = i + 1;  
end
x_axis = x_axis * step_size;


Btrad = alltrad;
Bsocial = allsocial;
Bind = allindividual;
Btot = total;
Bret = indret;
POPB = inddists;

inddists = (.6*POPB + .4*POPA) / 2;

figure;
hold on;
plot(x, inddists(12*30*10, :)/trapz(x, inddists(12*30*10, :)), 'r', "DisplayName", "12month Distribution");
plot(x, inddists(end,:) / trapz(x, inddists(end, :)), 'b', "DisplayName", "12month + 1day Distribution");
xlabel('Ideological scale');
ylabel('Normalized proportion of ideology distribution');
title("Ideological Distribution of individual")
legend show;


indret = (.6*Aret + .4*Bret)  ;
alltrad = (.6*Atrad + .4*Btrad);
allsocial = (.6*Asocial + .4*Bsocial);
allindividual = (.6*Aind + .4*Bind) ;

figure;
hold on;
plot(x_axis, alltrad(1:length(x_axis)), 'r', 'DisplayName', 'Traditional');
plot(x_axis, allsocial(1:length(x_axis)), 'b', 'DisplayName', 'Social');
plot(x_axis, allindividual(1:length(x_axis)), 'g', 'DisplayName', 'Individual');
plot(x_axis, total(1:length(x_axis)), 'DisplayName', 'Total');
plot(x_axis, indret(1:length(x_axis)), 'DisplayName', 'retained');

% Adding titles and labels
title('INF Over Time');
xlabel('Time (days)');
ylabel('Proportion of INF expected');
legend show;
grid on;

avg1 = sum(x .*  inddists(12*30*10, :)) / sum( inddists(12*30*10, :));

avg2 =  sum(x .*  inddists((12+1)*30*10, :)) / sum( inddists((12+1)*30*10, :));

disp(avg1);
disp(avg2)
