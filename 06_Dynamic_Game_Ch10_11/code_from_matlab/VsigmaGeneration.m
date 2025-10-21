function [W1out, W2out] = VsigmaGeneration(CCP1,CCP2, EstimatedTransition, ...
    EVrandom, UNIrandom, InitialState, NumSimMarkets, NumSimulations, NumSimPeriods, beta)

% Preparation
W1 = zeros(6,NumSimMarkets,NumSimulations);
W2 = zeros(6,NumSimMarkets,NumSimulations);

for mrkt = 1:NumSimMarkets

    %     % Extracting the data
    %     State = InitialState(mrkt,2);

    % The threshold values should not be varied by simulations
    % and these are 8*1 vectors
    ThresholdValue1 = log(CCP1) - log(ones(8,1)-CCP1);
    ThresholdValue2 = log(CCP2) - log(ones(8,1)-CCP2);

    for sim = 1:NumSimulations

        % Extracting the data
        State = InitialState(mrkt,2);

        for t=1:NumSimPeriods

            %             disp("period")
            %             disp(t)
            %             disp("state")
            %             disp(State)

            DiffEpsilon1inc = - EVrandom(mrkt,t,1,sim,State,2) + EVrandom(mrkt,t,1,sim,State,3);
            DiffEpsilon1dec = - EVrandom(mrkt,t,1,sim,State,2) + EVrandom(mrkt,t,1,sim,State,1);

            DiffEpsilon2inc = - EVrandom(mrkt,t,2,sim,State,2) + EVrandom(mrkt,t,2,sim,State,3);
            DiffEpsilon2dec = - EVrandom(mrkt,t,2,sim,State,2) + EVrandom(mrkt,t,2,sim,State,1);

            DiffFirm1inc = ThresholdValue1(State) - DiffEpsilon1inc;
            DiffFirm1dec = ThresholdValue1(State) - DiffEpsilon1dec;

            DiffFirm2inc = ThresholdValue2(State) - DiffEpsilon2inc;
            DiffFirm2dec = ThresholdValue2(State) - DiffEpsilon2dec;

            % a1 = 0 means stay, otherwise (a1 =1 ) means increasing or decreasing
            if (State == 1)||(State == 2)||(State == 5)||(State == 6)
                % n1 = 0
                a1 = 0*(DiffFirm1inc > 0) + 1*(DiffFirm1inc < 0);
                e1 = (1-a1)*EVrandom(mrkt,t,1,sim,State,2) + a1*EVrandom(mrkt,t,1,sim,State,3);
            else
                % n1 = 1
                a1 = 0*(DiffFirm1dec > 0) + 1*(DiffFirm1dec < 0);
                e1 = (1-a1)*EVrandom(mrkt,t,1,sim,State,2) + a1*EVrandom(mrkt,t,1,sim,State,1);
            end

            % Firm 2's action
            if (State == 1)||(State == 3)||(State == 5)||(State == 7)
                % n2 = 0
                a2 = 0*(DiffFirm2inc > 0) + 1*(DiffFirm2inc < 0);
                e2 = (1-a2)*EVrandom(mrkt,t,2,sim,State,2) + a2*EVrandom(mrkt,t,2,sim,State,3);
                %                 e2 = (1-a2)*EVrandom(mrkt,t,1,sim,State,2) + a2*EVrandom(mrkt,t,1,sim,State,3);

            else
                % n1 = 1
                a2 = 0*(DiffFirm2dec > 0) + 1*(DiffFirm2dec < 0);
                e2 = (1-a2)*EVrandom(mrkt,t,2,sim,State,2) + a2*EVrandom(mrkt,t,2,sim,State,1);
                %                 e2 = (1-a2)*EVrandom(mrkt,t,1,sim,State,2) + a2*EVrandom(mrkt,t,1,sim,State,1);

            end

            % Finally, decide the next state based on the current
            % state, a1, a2, and transition matrix
            if State == 1

                W1Seed = [ 0 0 0 0 a1 e1 ];
                W2Seed = [ 0 0 0 0 a2 e2 ];

                %                 if (a1==0)&&(a2==0)
                %                     NextState = 1;
                %                 elseif (a1==0)&&(a2==1)
                %                     NextState = 2;
                %                 elseif (a1==1)&&(a2==0)
                %                     NextState = 3;
                %                 else
                %                     NextState =4;
                %                 end
                NextState = 2*a1+a2+1;

                TransitionSeed = (UNIrandom(mrkt,t,sim) > EstimatedTransition(1,1));

                %NextState = (1+TransitionSeed)*NextState;
                NextState = 4*TransitionSeed + NextState;

            elseif State == 2

                W1Seed = [ 0 0 0 0 a1 e1 ];
                W2Seed = [ 1 0 1 a2 0 e2 ];

                %                 if (a1==0)&&(a2==0)
                %                     NextState = 2;
                %                 elseif (a1==0)&&(a2==1)
                %                     NextState = 1;
                %                 elseif (a1==1)&&(a2==0)
                %                     NextState = 4;
                %                 else
                %                     NextState = 3;
                %                 end

                NextState = 2+2*a1-a2;

                TransitionSeed = (UNIrandom(mrkt,t,sim) > EstimatedTransition(1,1));
                %NextState = (1+TransitionSeed)*NextState;
                NextState = 4*TransitionSeed + NextState;

            elseif State == 3

                W1Seed = [ 1 0 1 a1 0 e1 ];
                W2Seed = [ 0 0 0 0 a2 e2 ];

                %                 if (a1==0)&&(a2==0)
                %                     NextState = 3;
                %                 elseif (a1==0)&&(a2==1)
                %                     NextState = 4;
                %                 elseif (a1==1)&&(a2==0)
                %                     NextState = 1;
                %                 else
                %                     NextState = 2;
                %                 end

                NextState = 3-2*a1+a2;

                TransitionSeed = (UNIrandom(mrkt,t,sim) > EstimatedTransition(1,1));
                %NextState = (1+TransitionSeed)*NextState;
                NextState = 4*TransitionSeed + NextState;

            elseif State == 4

                W1Seed = [ 1 1 1 a1 0 e1 ];
                W2Seed = [ 1 1 1 a2 0 e2 ];

                %                 if (a1==0)&&(a2==0)
                %                     NextState = 4;
                %                 elseif (a1==0)&&(a2==1)
                %                     NextState = 3;
                %                 elseif (a1==1)&&(a2==0)
                %                     NextState = 2;
                %                 else
                %                     NextState = 1;
                %                 end

                NextState = 4-2*a1-a2;

                TransitionSeed = (UNIrandom(mrkt,t,sim) > EstimatedTransition(1,1));
                %NextState = (1+TransitionSeed)*NextState;
                NextState = 4*TransitionSeed + NextState;

            elseif State == 5

                W1Seed = [ 0 0 0 0 a1 e1 ];
                W2Seed = [ 0 0 0 0 a2 e2 ];

                %                 if (a1==0)&&(a2==0)
                %                     NextState = 1;
                %                 elseif (a1==0)&&(a2==1)
                %                     NextState = 2;
                %                 elseif (a1==1)&&(a2==0)
                %                     NextState = 3;
                %                 else
                %                     NextState =4;
                %                 end

                NextState = 1+2*a1+a2;

                TransitionSeed = (UNIrandom(mrkt,t,sim) < EstimatedTransition(2,2));
                %NextState = (1+TransitionSeed)*NextState;
                NextState = 4*TransitionSeed + NextState;

            elseif State == 6

                W1Seed = [ 0 0 0 0 a1 e1 ];
                W2Seed = [ 1 0 0 a2 0 e2 ];

                %                 if (a1==0)&&(a2==0)
                %                     NextState = 2;
                %                 elseif (a1==0)&&(a2==1)
                %                     NextState = 1;
                %                 elseif (a1==1)&&(a2==0)
                %                     NextState = 4;
                %                 else
                %                     NextState = 3;
                %                 end

                NextState = 2+2*a1-a2;

                TransitionSeed = (UNIrandom(mrkt,t,sim) < EstimatedTransition(2,2));
                %NextState = (1+TransitionSeed)*NextState;
                NextState = 4*TransitionSeed + NextState;


            elseif State == 7

                W1Seed = [ 1 0 0 a1 0 e1 ];
                W2Seed = [ 0 0 0 0 a2 e2 ];

                %                 if (a1==0)&&(a2==0)
                %                     NextState = 3;
                %                 elseif (a1==0)&&(a2==1)
                %                     NextState = 4;
                %                 elseif (a1==1)&&(a2==0)
                %                     NextState = 1;
                %                 else
                %                     NextState = 2;
                %                 end

                NextState = 3-2*a1+a2;

                TransitionSeed = (UNIrandom(mrkt,t,sim) < EstimatedTransition(2,2));
                %NextState = (1+TransitionSeed)*NextState;
                NextState = 4*TransitionSeed + NextState;

            elseif State == 8

                W1Seed = [ 1 1 0 a1 0 e1 ];
                W2Seed = [ 1 1 0 a2 0 e2 ];

                %                 if (a1==0)&&(a2==0)
                %                     NextState = 4;
                %                 elseif (a1==0)&&(a2==1)
                %                     NextState = 3;
                %                 elseif (a1==1)&&(a2==0)
                %                     NextState = 2;
                %                 else
                %                     NextState = 1;
                %                 end

                NextState = 4-2*a1-a2;

                TransitionSeed = (UNIrandom(mrkt,t,sim) < EstimatedTransition(2,2));
                %NextState = (1+TransitionSeed)*NextState;
                NextState = 4*TransitionSeed + NextState;

            else
                momomomo

            end

            %CheckStateDecision(NumSimPeriods*(sim-1)+t,1) = sim;
            %CheckStateDecision(NumSimPeriods*(sim-1)+t,2) = t;
            %CheckStateDecision(NumSimPeriods*(sim-1)+t,3) = State;
            %CheckStateDecision(NumSimPeriods*(sim-1)+t,4) = a1;
            %CheckStateDecision(NumSimPeriods*(sim-1)+t,5) = a2;

            W1(:,mrkt,sim) = W1(:,mrkt,sim) + beta^(t-1)*W1Seed';
            W2(:,mrkt,sim) = W2(:,mrkt,sim) + beta^(t-1)*W2Seed';

            %             [sim,t,a1,a2]

            State = NextState;
        end

    end
end

W1out = mean(W1(:,:,1:NumSimulations),3);
W2out = mean(W2(:,:,1:NumSimulations),3);
