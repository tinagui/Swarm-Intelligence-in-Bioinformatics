function [msp] = MaMove(fn,mn,fsp,msp,femass,mamass,d,lb,ub,pm)
%MAMOVE Summary of this function goes here
%   Detailed explanation goes here
    %Preliminaries
    dt=zeros(1,mn);
    % Scale for distance
    scale=(-lb(1)+ub(1));%/2;
    [Indb,~] = find(mamass>=median(mamass));  % Male spiders above median
    for i=1:mn
        if ismember(i,Indb)     % Spider above the median
            % Start looking for a female with stronger vibration
            for j=1:fn
                if femass(j)>mamass(i)
                    % Calculate the distance
                    dt(j)=norm(msp(i,:)-fsp(j,:));
                else
                    dt(j)=0;
                end
            end
            % Choose the shortest distance
            [~,Ind,val] = find(dt);   % Choose where the distance in non zero
            [~,Imin] = min(val);      % Get the shortest distance
            Ish = Ind(Imin);
            % Update moves
            if isempty(val)
                Vib=0;
                spaux=zeros(1,d);
            else
                dt=dt./scale;
                Vib = 2*femass(Ish)*exp(-(rand*dt(Ish).^2));
                spaux=fsp(Ish,:);
            end
            delta = 2*rand(1,d)-.5;
            tmpf = 2*pm.*(rand(1,d)-0.5);
            msp(i,:) = msp(i,:)+Vib*(spaux-msp(i,:)).*delta+tmpf;
        else % de aqui para abajo falta
            %% Spider below median, go to weigthed mean
            % Generate the weighted mean
            spdpos = [fsp' msp']';
            spdwei = [femass' mamass']';
            weigth = repmat(spdwei,1,d);
            dim = find(size(spdpos)~=1,1);
            wmean = sum(weigth.*spdpos,dim)./sum(weigth,dim);
            %% Move
            delta = 2*rand(1,d)-.5;
            tmpf = 2*pm.*(rand(1,d)-0.5);
            msp(i,:) = msp(i,:)+(wmean-msp(i,:)).*delta+tmpf;
        end
    end
    % Check limits
        for i=1:d
            for j=1:mn
                if msp(j,i)<lb(i), msp(j,i)=lb(i)+(ub(i)-lb(i)).*rand(1,1); end
                if msp(j,i)==lb(i), msp(j,i)=lb(i); end

                if msp(j,i)>ub(i), msp(j,i)=lb(i)+(ub(i)-lb(i)).*rand(1,1); end
                if msp(j,i)==ub(i), msp(j,i)=ub(i); end
            end
        end
end