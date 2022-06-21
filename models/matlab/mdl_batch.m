function [ e res ] = mdl_batch( X, batch )
%
%
%

% prepare batch data
if ~isa( batch, 'nominal' )
	batch = nominal( batch );
end
dbatch = double( batch );
bnames = getlevels( batch );
bn = numel( bnames );

% sizes
[ d N ] = size( X );

% paramas
bnin = 400;
nsamples = 400;
thin = 5;
[ ap am ] = deal( 0.1, 8 );
[ cp cm ] = deal( 0.1, 0 );
[ bs br ] = deal( 1.1, 0.001 );

% init
a = mean( X, 2 );
e = bsxfun( @minus, X, a );

c = zeros( d, bn );
for j=1:bn
	c(:,j) = mean( e(:,dbatch == j), 2 );
end

e = bsxfun( @minus, X - c(:,dbatch), a );
b = var( e, [], 2 );

% traces
ns = floor( nsamples/thin );
tr.a = zeros( d, ns );
tr.c = zeros( d, bn, ns );
tr.b = zeros( d, ns );

% loop
ss = 1;
for s=-bnin:nsamples
	% resample a
	e = X - c(:,dbatch);
	
	Sj = 1./( ap + N./b );
	mj = Sj.*( am*ap + sum( e, 2 )./b );
	a = mj + sqrt( Sj ).*randn( d, 1 );
	
	% resample c
	e = bsxfun( @minus, X, a );
	
	for j=1:bn
		ix = dbatch == j;
		Sj = 1./( cp + sum( ix )./b );
		mj = Sj.*( cm*cp + sum( e(:,ix), 2 )./b );
		c(:,j) = mj + sqrt( Sj ).*randn( d, 1 );
	end
		
	% resample d
	e = bsxfun( @minus, X - c(:,dbatch), a );
	b = ( br + 0.5*sum( e.^2, 2 ) )./randg( bs + 0.5*N, d, 1 );
	
	% verbosity
	if mod( s, 20 ) == 0
		fprintf( 'it: %d\n', s )
	end
	
	% traces
	if s > 0 && mod( s, thin ) == 0
		tr.a(:,ss) = a;
		tr.c(:,:,ss) = c;
		tr.b(:,ss) = b;
		ss = ss + 1;
	end
end

res.tr = tr;
res.a = median( tr.a, 2 );
res.c = median( tr.c, 3 );
res.b = median( tr.b, 2 );

e = bsxfun( @minus, X - res.c(:,dbatch), res.a );
% e = bsxfun( @rdivide, bsxfun( @minus, X - res.c(:,dbatch), res.a ), sqrt( res.b ) );

end
