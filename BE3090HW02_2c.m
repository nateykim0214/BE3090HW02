%%%%%%%%%%%%%%%%% DISTRIBUTIONS %%%%%%%%%%%%%%%%%%%%%%%%
x = 0:.1:15;
basic = x;
traddist = normpdf(x, 5, 1);

s1 = normpdf(x, 3, .5);
s2 = normpdf(x, 7, .5);
scombined = s1+s2;
socialdist = scombined / trapz(x,scombined);
%3aa united individuals
inddist = normpdf(x, 5.5, 1);

%3ab split individuals
%p1 = normpdf(x, 6, 1);
%p2 = normpdf(x, 4, 1);
%combined = p1+p2;
%inddist = combined / trapz(x, combined);

%%%%%%%%%%%%%%%% END DISTRIBUTIONS %%%%%%%%%%%%%%%%%%%%%

misinfdist = normpdf(x, 9, 0.0000001);
BETA = 0.00000020205;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
months = 1;
days = months * 30;
t = (days+1)*10; %(days * 10) 
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


injecttime = 5;
injectamount = .05;

tradists = zeros(total_length, length(basic));
tradists(1, :) = traddist; 

socialdists = zeros(total_length, length(basic)); 
socialdists(1, :) = socialdist; 

inddists = zeros(total_length, length(basic));  
inddists(1, :) = inddist; 

totaldists = zeros(total_length, length(basic));  
totaldists(1, :) = traddist + socialdist + inddist;


x_axis = [0];  % Will store the days incremented by whole numbers

alltrad = [0.00001];  % Initialize traditional array
allsocial = [0.00001];  % Initialize social array
allindividual = [0.00001];  % Initialize individual array

mistrad = [0.000000001];
missocial = [0.00000001];
misindividual = [0.000000001];

indret = [0.00001];

current_x = 0;  % Current position in time
i = 1;  % Index for the arrays
total = [0];

kreti = .05003;
krettrad = kreti/3;
kretsoc = kreti/6;

tol = .005;  % Define a small tolerance

%%%%%%%%%%%%%%%%%%%%%%% GENERATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% While loop for calculations
while current_x < t
    % Convert current time to days for x-axis
    if mod(i-1, (1/step_size)) == 0  % Every 100 steps, increment the day
        x_axis(end + 1) = (i-1) / (1/step_size);  % Increment day by 1
    end


    if i*step_size  == injecttime+.2
        disp("Inject")
        alltrad(i) = alltrad(i) + injectamount/2;
        allsocial(i) = allsocial(i) + injectamount/2;
        mistrad(i) = injectamount/2;
        missocial(i) = injectamount/2;
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

    newtraddist = tradists(i, :);
    newsocialdist = socialdists(i, :);
    newinddist = inddists(i, :);

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

        if abs(newmistrad) > tol
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
        
        if abs(newmissocial) > tol
            newmistrad = newmistrad + changeFromFlow*allsocial(i);
            newmissocial = newmissocial - changeFromFlow*allsocial(i);
        end

        newtraddist(1, :) = newtraddist(1, :) + changeFromFlow*socialdists(i, :) + changeTfromGrowth * tradists(i, :);
        newsocialdist(1, :) =newsocialdist(1, :) - changeFromFlow*socialdists(i, :) + changeSfromGrowth * socialdists(i, :);

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
        
        if abs(newmissocial) > tol
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
        
        if abs(newmistrad) > tol
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
    
   
    % Store the updated values in the arrays
    alltrad(i + 1) = new_trad;
    allsocial(i + 1) = new_social;
    allindividual(i + 1) = new_individual;
    total(i+1) = new_trad + new_social + new_individual;
    
    indret(i+1) = indret(i)+ (new_trad * krettrad + new_social * kretsoc + new_individual * kreti)*step_size;

    tradists(i+1, :) = newtraddist(1, :);
    socialdists(i+1, :) = newsocialdist(1, :);
    inddists(i+1, :) = newinddist(1, :);
    totaldists(i+1, :) = newtraddist(1, :) + newsocialdist(1, :) + newinddist(1, :);

    mistrad(i+1) = newmistrad;
    missocial(i+1) =  newmissocial;
    misindividual(i+1) = newmisindividual;
    % Update current time
    current_x = current_x + step_size;  % Increment time
    i = i + 1;  % Increment index
end

x_axis = x_axis * step_size;


% Plotting alltrad, allsocial, and allindividual vs x_axis
figure;
hold on;
plot(x_axis, alltrad(1:length(x_axis)), 'r', 'DisplayName', 'Traditional');
plot(x_axis, allsocial(1:length(x_axis)), 'b', 'DisplayName', 'Social');
plot(x_axis, allindividual(1:length(x_axis)), 'g', 'DisplayName', 'Individual');
plot(x_axis, total(1:length(x_axis)), 'DisplayName', 'Total');
plot(x_axis, mistrad(1:length(x_axis)), ':', 'LineWidth', 1.5, 'DisplayName', "trad mis");
plot(x_axis, missocial(1:length(x_axis)), ':', 'LineWidth', 1.5,'DisplayName', "social mis");
plot(x_axis, misindividual(1:length(x_axis)), ':', 'LineWidth', 1.5,'DisplayName', "individual mis");
hold off;

% Adding titles and labels
title('INF Over Time (once every 7 days)');
xlabel('Time (days)');
ylabel('Proportion of INF expected');
legend show;
grid on;