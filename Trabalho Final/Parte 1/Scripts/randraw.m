
function varargout = randraw(distribName, distribParams, varargin)
		funcName = mfilename;
		if nargin == 0
				help(funcName);
				return;
		elseif nargin == 1 runMode = 'distribHelp';
		elseif nargin == 2
				runMode = 'genRun';
				sampleSize = [1 1];
		else
				runMode = 'genRun';
				sampleSize = [varargin{1:end}];
		end

		distribNameInner = lower( distribName( ~isspace( distribName ) ) );
		if strcmp(runMode, 'distribHelp')
				fid = fopen( [ funcName '.m' ], 'r' );
				printHelpFlag = 0;
				while 1
						tline = fgetl( fid );
						if ~ischar( tline )
								fprintf( '\n Unknown distribution name ''%s''.\n', distribName );
								break;
						end

						if ~isempty( strfind( tline, [ 'END ', distribNameInner,' HELP' ] ) )
								printHelpFlag = 0;
								break;
						end

						if printHelpFlag
								startPosition = strfind( tline, ' % ' ) + 3;
								printLine = tline( startPosition : end );
								if ~strcmp( funcName, 'randraw' )
										indxs = strfind( printLine, 'randraw' );
										while ~isempty( indxs )
												headLine = printLine( 1:indxs(1)-1 );
												tailLine = printLine( indxs(1)+7:end );
												printLine = [ headLine, funcName, tailLine ];
												indxs = strfind( printLine, 'randraw' );
										end
								end, pause(0.02);
								fprintf( '\n%s', printLine );
						end

						if ~isempty( strfind( tline, [ 'START ', distribNameInner,' HELP' ] ) ) printHelpFlag = 1; end
				end, fprintf( '\n\n' ); fclose( fid );

				if nargout > 0 varargout{1} = []; end
				return;
		end

		if length(sampleSize) == 1 sampleSize = [ sampleSize, 1 ]; end
		if strcmp(runMode, 'genRun')
				runExample = 0;
				plotFlag = 0;

				dbclear if warning;
				out = [];
				 if prod(sampleSize) > 0
						switch lower( distribNameInner )
								case {'gamma'}
										checkParamsNum(funcName, 'Gamma', 'gamma', distribParams, [1, 3]);
										if numel(distribParams)==3
												m  = distribParams(1);
												b  = distribParams(2);
												a  = distribParams(3);
												validateParam(funcName, 'Gamma', 'gamma', '[m, b, a]', 'a', a, {'> 0'});
												validateParam(funcName, 'Gamma', 'gamma', '[m, b, a]', 'b', b, {'> 0'});
										else
												m = 0;
												b = 1;
												a  = distribParams(1);
												validateParam(funcName, 'Gamma', 'gamma', '[m, b, a]', 'a', a, {'> 0'});
										end

										if a < 1
												out = m + b*(feval(funcName, 'gamma', 1+a, sampleSize)).*(rand(sampleSize).^(1/a));
										else
												d = a - 1/3;
												c = 1/sqrt(9*d);
												x = randn( sampleSize );
												v = 1+c*x;

												indxs = find(v <= 0);
												while ~isempty(indxs)
														indxsSize = size( indxs );
														xNew = randn( indxsSize );
														vNew = a+c*xNew;

														l = (vNew > 0);
														v( indxs( l ) ) = vNew(l);
														x( indxs( l ) ) = xNew(l);
														indxs = indxs( ~l );
												end

												u = rand( sampleSize );
												v = v.^3;
												x2 = x.^2;
												out = d*v;

												indxs = find( (u>=1-0.0331*x2.^2) & (log(u)>=0.5*x2+d*(1-v+log(v))) );
												while ~isempty(indxs)
														indxsSize = size( indxs );

														x = randn( indxsSize );
														v = 1+c*x;
														indxs1 = find(v <= 0);
														while ~isempty(indxs1)
																indxsSize1 = size( indxs1 );
																xNew = randn( indxsSize1 );
																vNew = a+c*xNew;

																l1 = (vNew > 0);
																v( indxs1(l1) ) = vNew(l1);
																x( indxs1(l1) ) = xNew(l1);
																indxs1 = indxs1( ~l1 );
														end

														u = rand( indxsSize );
														v = v .* v .* v;
														x2 = x.*x;

														l = (u<1-0.0331*x2.*x2) | (log(u)<0.5*x2+d*(1-v+log(v)));
														out( indxs( l ) ) = d*v(l);
														indxs = indxs( ~l );
												end, out = m + b*out;
										end

								case {'nakagamim', 'nakagami-m', 'nakagami'}
										checkParamsNum(funcName, 'Nakagami', 'nakagamim', distribParams, [1 2]);
										if numel(distribParams) == 2
												m = distribParams(1);
												omega = distribParams(2);
												validateParam(funcName, 'Nakagami', 'nakagamim', '[m, omega]', 'm', m, {'> 0'});
												validateParam(funcName, 'Nakagami', 'nakagamim', '[m, omega]', 'omega', omega, {'> 0'});
										else
												m = distribParams(1);
												omega = 1;
												validateParam(funcName, 'Nakagami', 'nakagamim', '[m, omega]', 'm', m, {'> 0'});
										end, out = sqrt( randraw('gamma', [0, omega/m, m], sampleSize) );

								otherwise fprintf('\n RANDRAW: Unknown distribution name: %s \n', distribName);
						end
				end, varargout{1} = out;
				return;
		end, return;


function checkParamsNum(funcName, distribName, runDistribName, distribParams, correctNum)
		if ~any( numel(distribParams) == correctNum )
				error('%s Variates Generation:\n %s%s%s%s%s', ...
						distribName, 'Wrong numebr of parameters (run ', ...
						funcName, ...
						'(''', ...
						runDistribName, ...
						''') for help) ');
		end, return;

function validateParam(funcName, distribName, runDistribName, distribParamsName, paramName, param, conditionStr)
		condLogical = 1;
		eqCondStr = [];
		for nn = 1:length(conditionStr)
				if nn==1 eqCondStr = [eqCondStr conditionStr{nn}];
				else eqCondStr = [eqCondStr ' and ' conditionStr{nn}]; end

				eqCond = conditionStr{nn}(1:2);
				eqCond = eqCond(~isspace(eqCond));
				switch eqCond
						case{'<'} condLogical = condLogical & (param<str2num(conditionStr{nn}(3:end)));
						case{'<='} condLogical = condLogical & (param<=str2num(conditionStr{nn}(3:end)));
						case{'>'} condLogical = condLogical & (param>str2num(conditionStr{nn}(3:end)));
						case{'>='} condLogical = condLogical & (param>=str2num(conditionStr{nn}(3:end)));
						case{'~='} condLogical = condLogical & (param~=str2num(conditionStr{nn}(3:end)));
						case{'=='}
								if strcmp(conditionStr{nn}(3:end),'integer') condLogical = condLogical & (param==floor(param));
								else condLogical = condLogical & (param==str2num(conditionStr{nn}(3:end))); end
				end
		end

		if ~condLogical
				error('%s Variates Generation: %s(''%s'',%s, SampleSize);\n Parameter %s should be %s\n (run %s(''%s'') for help)', distribName, funcName, runDistribName, distribParamsName, ...
				paramName, eqCondStr, funcName, runDistribName); end, return;

