classdef ShimOptAcdc < ShimOpt
    
    methods
        
        function Shim = ShimOptAcdc( Params )
            if nargin < 1 || isempty( Params )
                Params.dummy = [] ;
            end
        end
        
        function Shim = optimizeshimcurrents( Shim, Params )
            Specs = ShimSpecs( Params );
            
            if nargin < 1
                error('Function requires at least 1 argument of type ShimOpt')
            elseif nargin == 1
                Params.dummy = [];
            end
            
            
            if ~myisfield(Params, 'maxCurrentPerChannel') || isempty( Params.maxCurrentPerChannel )
                Params.maxCurrentPerChannel = Specs.Amp.maxCurrentPerChannel ;
            end
            
            A = Shim.getshimoperator ;
            M = Shim.gettruncationoperator ;
            
            b = M*(-Shim.Field.img(:)) ;
            
            % -------
            % Least-squares solution via conjugate gradients
            Shim.Model.currents = cgls( A'*M'*M*A, ... % least squares operator
                A'*M'*b, ... % effective solution vector
                zeros( [Shim.Params.nActiveChannels 1] ), ... % initial model (currents) guess
                Shim.Params.CG ) ;
            % -------
            
            isCurrentSolutionOk = false ;
            
            if abs(Shim.Model.currents) <= Params.maxCurrentPerChannel
                isCurrentSolutionOk = true ;
                
            end
            
            if ~isCurrentSolutionOk
                
                Options = optimset(...
                    'DerivativeCheck','off',...
                    'GradObj','on',...
                    'Display', 'off',... %'iter-detailed',...
                    'MaxFunEvals',36000,...
                    'TolX',1e-11,...
                    'TolCon',1E-8);
                
                A = M*A ;
                
                
                tic
                [Shim.Model.currents] = fmincon( ...
                    @shimcost,...
                    ones(Specs.Amp.nActiveChannels,1),...
                    [],...
                    [],...
                    [],...
                    [],...
                    -Params.maxCurrentPerChannel*ones(Specs.Amp.nActiveChannels,1),...
                    Params.maxCurrentPerChannel*ones(Specs.Amp.nActiveChannels,1),...
                    Options);
                toc
            end
            
            function [f, df] = shimcost( currents )
                
                y=A*currents - b;
                f=y'*y;
                df=2*A'*y;
            end
            
            Shim = Shim.setforwardmodelfield ;
            
        end
    end
end