%%%%%%%%%%%%%%%%% DISTRIBUTIONS %%%%%%%%%%%%%%%%%%%%%%%%
x = 0:.5:10;
basic = x;
traddist = normpdf(x, 5, 1);

socialdist = (normpdf(x, 3, .5) + normpdf(x, 7, .5)) ;

inddist = oldinddist;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISINF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


misinfdist = normpdf(x, 1, 0.001) ;
%1,60.3; 2,60.4; 3, 0.45; 60.48;
BETA = 60.48;
BETAFLAG = false;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
months = 12;
days = months * 30;
t = days*10; %(days * 10) 
step_size = 0.1;  % Step size
total_length = round(t / step_size);

pts = 0.576;    
pst = 2.302;    
pti = 0.329;  
psi = 0.329;    

klostrad = .657;
tradperc = 1;
klostsocial = 1.53;
socialperc = 1;
klostind = .153;
indperc = 1;

kgsocial = .45;
kgind = .1;
kgtrad = .45;

socinsoc = [oldsocinsoc];
tradintrad = [oldtradintrad];
indinind = [oldindinind];


tradists = zeros(total_length, length(basic));
tradists(1, :) = oldtraddist; 

socialdists = zeros(total_length, length(basic)); 
socialdists(1, :) = socialdist; 

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

bigX = 1.16;
kreti =.05003;
krettrad = kreti/bigX;
kretsoc = kreti/6;
kretmis = kreti/2;


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



    newmistrad = mistrad(i);
    newmissocial = missocial(i);
    newmisindividual = misindividual(i);
    new_mis = newmistrad + newmissocial + newmisindividual;

    newmisindividual = newmisindividual + (- klostmis * newmisindividual) * step_size;
    new_mis = new_mis + (- klostmis * new_mis) * step_size;
    newmissocial = newmissocial + (- klostmis * newmissocial) * step_size;
    newmistrad = newmistrad + (- klostmis * newmistrad) * step_size;

    newtraddist = tradists(i, :);
    newsocialdist = socialdists(i, :);
    newinddist = inddists(i, :);


    newsocinsoc = socinsoc(i);
    newtradintrad = tradintrad(i);
    newindinind = indinind(i);
    new_social = new_social + (kgsocial - klostsocial * newsocinsoc)* socialperc *step_size;
        newsocinsoc = newsocinsoc + (kgsocial - klostsocial * newsocinsoc)* socialperc *step_size;
    new_trad = new_trad + (kgtrad  - klostrad * newtradintrad) * tradperc*step_size ;
        newtradintrad = newtradintrad + (kgtrad  - klostrad * newtradintrad) * tradperc*step_size ;
    new_individual = new_individual + (kgind  - klostind * newindinind)* indperc*step_size ;
        newindinind = newindinind +  (kgind  - klostind * newindinind)* indperc*step_size ;

    changeI = (kgind  - klostind * newindinind)* indperc  *step_size/ allindividual(i);
    changeT = (kgtrad  - klostrad * newtradintrad)* tradperc *step_size / alltrad(i);
    changeS = (kgsocial  - klostsocial * newsocinsoc)* socialperc  * step_size/ allsocial(i);
    
    newtraddist = newtraddist + (changeT)*traddist;
    newinddist = newinddist + (changeI)*newinddist;
    newsocialdist = newsocialdist + (changeS)*socialdist;
   

    %Misinfinjection
    if mod((i*step_size), 12*30) == 0 && BETAFLAG == true % @ 12 months
        misinfdist= misinfdist * BETA * .01 * trapz(x, totaldists(i, :));

        newtraddist = newtraddist + misinfdist/2 ;
        newsocialdist = newsocialdist  + misinfdist/2;

        new_trad = new_trad + BETA  * .01* total(i)/2;
        new_social = new_social + BETA  * .01* total(i)/2;

        newmistrad = mistrad(i) +  BETA * .01 * total(i)/2;
        newmissocial = missocial(i) + BETA * .01 * total(i)/2;
        

        
        BETAFLAG = false;
    end
 

    % Update for traditional vs social
    if diffTS > 0
        changeFromFlow = (pts * diffTS)  * step_size / alltrad(i);

        newtraddist = newtraddist - changeFromFlow * tradists(i,:);
        newsocialdist = newsocialdist + changeFromFlow * tradists(i, :);
    
        new_trad = new_trad - (pts * diffTS) * step_size;
            newtradintrad = newtradintrad - (pts*diffTS)*step_size;
        new_social = new_social + (pts * diffTS) * step_size;
        
        newmistrad = newmistrad - changeFromFlow * newmistrad;
        newmissocial = newmissocial + changeFromFlow * newmistrad;

    elseif diffTS < 0
        changeFromFlow = (pst * -diffTS) * step_size  / allsocial(i);

        newtraddist = newtraddist + changeFromFlow*socialdists(i, :);
        newsocialdist = newsocialdist - changeFromFlow*socialdists(i, :);

        new_trad = new_trad + (pst * -diffTS) * step_size;
        new_social = new_social - (pst*-diffTS) * step_size;

        newmissocial = newmissocial - changeFromFlow*newmissocial;
        newmistrad = newmistrad + changeFromFlow*newmistrad;
    end
    
    if diffIS < 0
        changeFromFlow = (psi * -diffIS)* step_size  / allsocial(i);

        newinddist = newinddist + changeFromFlow*socialdists(i, :);
        newsocialdist = newsocialdist - changeFromFlow*socialdists(i, :);

        new_individual = new_individual + (psi * -diffIS) * step_size;
        new_social = new_social - (psi * -diffIS) * step_size;

        newmisindividual = newmisindividual + changeFromFlow * newmissocial;
        newmissocial = newmissocial - changeFromFlow * newmissocial;
    end

    if diffTI > 0
        changeFromFlow = (pti * diffTI) * step_size  / alltrad(i);

        newinddist = newinddist + changeFromFlow*tradists(i,:) ;
        newtraddist = newtraddist - changeFromFlow*tradists(i,:) ;
       
        new_trad = new_trad - (pti * diffTI)* step_size;
        new_individual = new_individual + (pti * diffTI)* step_size ;

        newmisindividual = newmisindividual + changeFromFlow * newmistrad;
        newmistrad = newmistrad - changeFromFlow *newmistrad;
    end


    socinsoc(i+1) = newsocinsoc;
    tradintrad(i+1) = newtradintrad;
    indinind(i+1) = newindinind;

    alltrad(i + 1) = new_trad;
    allsocial(i + 1) = new_social;
    allindividual(i + 1) = new_individual;
    total(i+1) = new_trad + new_social + new_individual;

    mistrad(i+1) = newmistrad;
    missocial(i+1) =  newmissocial;
    misindividual(i+1) = newmisindividual;
    allmis(i+1) = newmistrad + newmissocial + newmisindividual;
    
    indret(i+1) = indret(i)+ (new_trad * krettrad + new_social * kretsoc + new_individual * kreti)*step_size;
    
    tradists(i+1, :) = newtraddist(1, :);
    socialdists(i+1, :) = newsocialdist(1, :);
    inddists(i+1, :) = newinddist(1, :);

    totaldists(i+1, :) = newtraddist(1, :) + newsocialdist(1, :) + newinddist(1, :);

    current_x = current_x + step_size; 
    i = i + 1;  
        % figure;
        % hold on;
        % plot(x, newtraddist(1,:))
        % plot(x, newsocialdist(1,:))
        % plot(x, totaldists(i,:))
        % legend show;
end


x_axis = x_axis * step_size;

% Plotting alltrad, allsocial, and allindividual vs x_axis
figure;
hold on;
plot(x_axis, alltrad(1:length(x_axis)), 'r', 'DisplayName', 'Traditional');
plot(x_axis, allsocial(1:length(x_axis)), 'b', 'DisplayName', 'Social');
plot(x_axis, allindividual(1:length(x_axis)), 'g', 'DisplayName', 'Individual');
plot(x_axis, total(1:length(x_axis)), 'DisplayName', 'Total');
plot(x_axis, indret(1:length(x_axis)), 'DisplayName', 'retained');
hold on;

title('INF Over Time');
xlabel('Time (days)');
ylabel('Proportion of INF expected');
legend show;
grid on;


figure;
hold on;
plot(x, inddists(1, :) / trapz(x, inddists(1, :)), 'b', "DisplayName", "end of a (12m)");
plot(x, inddists(6*30*10, :) / trapz(x, inddists(6*30*10, :)), 'b', "DisplayName", "a+6months (18m)");
plot(x, inddists(12*30*10, :) /trapz(x, inddists(12*30*10, :)) , 'b', "DisplayName", "a+12months (24m)");
xlabel('Ideological scale');
ylabel('Normalized proportion of ideology distribution');
title("Ideological Distribution of individual")
legend show;


oldtrad = alltrad(length(x_axis));
oldsocial = allsocial(length(x_axis));
oldindividual = allindividual(length(x_axis));
oldtotal = total(length(x_axis));
oldret = indret(length(x_axis));
disp(oldret)

oldtraddist = tradists(12*30*10, :); 
oldsocialdist = socialdists(12*30*10, :);
oldinddist = inddists(12*30*10, :);
oldtotaldist = totaldists(12*30*10, :);

oldsocinsoc = socinsoc(end);
oldtradintrad = tradintrad(end);
oldindinind = indinind(end);
