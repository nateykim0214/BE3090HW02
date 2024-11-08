%%%%%%%%%%%%%%%%% DISTRIBUTIONS %%%%%%%%%%%%%%%%%%%%%%%%
x = 0:.5:15;
basic = x;
traddist = normpdf(x, 5, 1);

socialdist = (normpdf(x, 3, .5) + normpdf(x, 7, .5)) / trapz(x, (normpdf(x, 3, .5) + normpdf(x, 7, .5)));

p1 = normpdf(x, 6, 1);
p2 = normpdf(x, 4, 1);
combined = p1+p2;
inddist = combined / trapz(x, combined);


%inddist = normpdf(x, 5.5, 1);
first = inddist;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISINF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


misinfdist = normpdf(x, 1, 0.000000000001) ;
misinfdist = misinfdist /  trapz(misinfdist);

BETA =  74;
BETAFLAG = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
months = 24;
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

tradists = zeros(total_length, length(basic));
tradists(1, :) = traddist; 

socialdists = zeros(total_length, length(basic)); 
socialdists(1, :) = socialdist; 


inddists = zeros(total_length, length(basic));  
inddists(1, :) = inddist; 

totaldists = zeros(total_length, length(basic));  
totaldists(1, :) = traddist + socialdist + inddist;

x_axis = [0];  % Will store the days incremented by whole numbers
alltrad = [0.001];  % Initialize traditional array
allsocial = [0.001];  % Initialize social array
allindividual = [0.001];  % Initialize individual array
allmis = [0.00000000000001];


mistrad = [0.0000000000001];
missocial = [0.00000000001];
misindividual = [0.000000000001];


indret = [0.001];
current_x = 0;  % Current position in time
i = 1;  % Index for the arrays
total = [.00000000000000000001];

bigX = 1.16;
kreti =.05003;
krettrad = kreti/bigX;
kretsoc = kreti/6;
kretmis = kreti/2;
C = 1000000;


kr = .1;
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

    sprimet= 6.94e-4  *step_size * .01 *(mistrad(i) / alltrad(i))^-1; 
    sprimes = 6.94e-4*step_size * .01 *(missocial(i) / allsocial(i))^-1; 
    sprimei = 6.94e-4 *step_size * .01 *(misindividual(i) / allindividual(i))^-1; 
        precogt = 1 - exp((-1/C) * (current_x / sprimet));
        precogs = 1 - exp((-1/C) * (current_x / sprimes));
        precogi = 1 - exp((-1/C) * (current_x / sprimei));  
            newmistrad = mistrad(i) - precogt*kr*mistrad(i) *step_size;
            newmissocial =  missocial(i) - precogs*kr*missocial(i)*step_size;
            newmisindividual = misindividual(i) - precogi*kr*misindividual(i)*step_size;

     
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

        newtraddist = newtraddist - changeFromFlow * (tradists(i,:));
        newsocialdist = newsocialdist + changeFromFlow * (tradists(i,:));
    
        new_trad = new_trad - (pts * diffTS) * step_size;
        new_social = new_social + (pts * diffTS) * step_size;
        
        newmistrad = newmistrad - changeFromFlow * newmistrad;
        newmissocial = newmissocial + changeFromFlow * newmistrad;
    elseif diffTS < 0
        changeFromFlow = (pst * -diffTS) * step_size  / allsocial(i);

        newtraddist = newtraddist + changeFromFlow*(socialdists(i,:)+0);
        newsocialdist = newsocialdist - changeFromFlow*(socialdists(i,:)+0);

        new_trad = new_trad - (pts * diffTS) * step_size;
        new_social = new_social + (pts * diffTS) * step_size;
        
        newmistrad = newmistrad - changeFromFlow * newmistrad;
        newmissocial = newmissocial + changeFromFlow * newmistrad;
    end
    if diffIS < 0
        changeFromFlow = (psi * -diffIS)* step_size  / allsocial(i);

        newinddist = newinddist + changeFromFlow*(socialdists(i,:)+0);
        newsocialdist = newsocialdist - changeFromFlow*(socialdists(i,:)+0);

        new_individual = new_individual + (psi * -diffIS) * step_size;
        new_social = new_social - (psi * -diffIS) * step_size;

        newmisindividual = newmisindividual + changeFromFlow * newmissocial;
        newmissocial = newmissocial - changeFromFlow * newmissocial;
    end
    if diffTI > 0
        changeFromFlow = (pti * diffTI) * step_size  / alltrad(i);

        newinddist = newinddist + changeFromFlow*(tradists(i,:)+0) ;
        newtraddist = newtraddist - changeFromFlow*(tradists(i,:)+0) ;
       
        new_trad = new_trad - (pti * diffTI)* step_size;
        new_individual = new_individual + (pti * diffTI)* step_size ;

        newmisindividual = newmisindividual + changeFromFlow * newmistrad;
        newmistrad = newmistrad - changeFromFlow *newmistrad;
    end

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
end
x_axis = x_axis * step_size;


% Plotting alltrad, allsocial, and allindividual vs x_axis
figure;
hold on;
plot(x_axis(3600:end), alltrad(3600:length(x_axis)), 'r', 'DisplayName', 'Traditional');
plot(x_axis(3600:end), allsocial(3600:length(x_axis)), 'b', 'DisplayName', 'Social');
plot(x_axis(3600:end), allindividual(3600:length(x_axis)), 'g', 'DisplayName', 'Individual');
plot(x_axis(3600:end), total(3600:length(x_axis)), 'DisplayName', 'Total');
%plot(x_axis(3600:end), indret(3600:length(x_axis)), 'DisplayName', 'retained');
plot(x_axis(3600:end), mistrad(3600:length(x_axis)), ':', 'LineWidth', 1.5, 'DisplayName', "trad mis");
plot(x_axis(3600:end), missocial(3600:length(x_axis)), ':', 'LineWidth', 1.5,'DisplayName', "social mis");
plot(x_axis(3600:end), misindividual(3600:length(x_axis)), ':', 'LineWidth', 1.5,'DisplayName', "individual mis");
plot(x_axis(3600:end), allmis(3600:length(x_axis)), ':', 'LineWidth', 1.5,'DisplayName', "all mis");

disp(length(x_axis))
for i = 3601:length(x_axis)
    if allmis(i)/total(i) <= .01
        disp(allmis(3601))
        disp(allmis(i))
        disp(i/10 - 360)
        break
    end
end

hold on;
% Adding titles and labels
title('INF Over Time');
xlabel('Time (days)');
ylabel('Proportion of INF expected');
legend show;
grid on;


figure;
hold on;
plot(x, inddists(end, :) / trapz(x, inddists(end, :)), 'r', "DisplayName", "Final Distribution");
plot(x, inddists((12+1)*30*10, :) /trapz(x, inddists((12+1)*30*10, :)) , 'g', "DisplayName", "Injection Distribution");
xlabel('Ideological scale');
ylabel('Normalized proportion of ideology distribution');
title("Ideological Distribution of individual")
legend show;


finalA = inddists;





