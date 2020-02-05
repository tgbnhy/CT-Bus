% call the functions
%%%%%%%%%%%% parameters
weights = [0.3, 0.5, 0.7]; % the weight to balance each dimension
ks = [10, 30, 50]; % number of new edges in a path
number_of_turnss = [1, 3, 5];% number of turns, candidates
seeding_numbers = [3000, 5000, 7000];% number of candidate edges as initial paths

%{
for weight = weights
    fairbus(weight, 30, 3, 5000);
end

for k = ks
    if k==30
        continue;
    end
    fairbus(0.5, k, 3, 5000);
end

for number_of_turns = number_of_turnss
    if number_of_turns == 3
        continue;
    end
    fairbus(0.5, 30, number_of_turns, 5000);
end

for seeding_number = seeding_numbers
    if seeding_number==5000
        continue;
    end
    fairbus(0.5, 30, 3, seeding_number);
end
%}
fairbus(0.5, 80, 3, 5000);