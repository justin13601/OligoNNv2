function exp_geo()
%
% Build the classification model presented in:
%
% Host gene expression classifiers diagnose acute respiratory illness etiology
% Sci Transl Med 2016
% PMID: 26791949
%
% Note: results may change slightly due to architecture/os differences, matlab
% and toolbox versions.
%

rng( 0 )

close all
clc

% get data from geo
rd = geoseriesread('GSE63990_series_matrix.txt');

% process key
batch = str2double( arrayfun( @(i) regexp( [ rd.Header.Samples.characteristics_ch1{:,i} ], '(?<=batch: )[01]', 'match' ), 1:size(rd.Header.Samples.characteristics_ch1, 2 ) ) )';
status = arrayfun( @(i) regexp( [ rd.Header.Samples.characteristics_ch1{:,i} ], '(?<=infection_status: )[\w \-]+', 'match' ), 1:size(rd.Header.Samples.characteristics_ch1, 2 ) )';
tmp = arrayfun( @(i) regexp( [ rd.Header.Samples.characteristics_ch1{:,i} ], '(?<=replicates: )\d', 'match' ), 1:size(rd.Header.Samples.characteristics_ch1, 2 ), 'uniformoutput', false );
idx = find( ~cellfun( @isempty, tmp ) );
repl = nan( size( tmp' ) );
repl(idx) = str2double( [ tmp{idx} ]' );

% batch correction
Xs = batch_filter( rd.Data, batch, repl );
% drop replicates from key
ix = batch == 0 & ~isnan( repl );
batch(ix) = []; status(ix) = [];

% show some numbers
fprintf( 'batches:\n' ), tabulate( batch )
fprintf( 'status:\n' ), tabulate( status )

% classifiers
[ ~, nzb ~, ypb ] = fit( Xs, nominal( status ), 'bacterial' );
[ ~, nzv ~, ypv ] = fit( Xs, nominal( status ), 'viral' );
[ ~, nzn ~, ypn ] = fit( Xs, nominal( status ), 'non-infectious illness' );

nz = [ nzb nzv nzn ];
pp = [ ypb ypv ypn ];

% display nz
fprintf( 'bacterial: %d viral: %d non-infectious illness: %d\n', sum( nz ) )

% compute confusion matrix
[ ~, ypm ] = max( pp, [], 2 );
cmt = confusionmat( double( reorderlevels( nominal( status ), { 'bacterial' 'viral' 'non-infectious illness' } ) ), ypm );
display( cmt )
display(  bsxfun( @rdivide, cmt, sum( cmt, 2 ) ) )

end

function [ oo nz id yp ] = fit( X, y, pos_class )
%
% Fit sparse logistic regression with leave-one-out cross-validation.
%
% Get cvglmnet from: http://web.stanford.edu/~hastie/glmnet_matlab/
%
% Returns: signature (nz), predisctions (yp), AUC, TPR and TNR.
%
% Note: this is not nested cross-validation, which takes approximately N
% (sample size) times longer to run, thus results may differ slighlty from
% those reported in the publication.
%

yb = y == pos_class;
oo = cvglmnet( X', yb, 'binomial', struct( 'nlambda', 50 ), 'class', size( X, 2 ), [], false, true, false );

[ auc thr tpr tnr ] = deal( zeros( size( oo.fit_preval, 2 ), 1 ) );
for n=1:size( oo.fit_preval, 2 )
	[ xx yy tt auc(n) ] = perfcurve( yb, oo.fit_preval(:,n), true );
	[ ~, iy ] = min( xx.^2 + ( yy - 1 ).^2 );
	thr(n) = tt(iy);
	tpr(n) = yy(iy);
	tnr(n) = 1 - xx(iy);
end

[ ~, id ] = max( auc );
[ xx yy tt auc_ ] = perfcurve( yb, oo.fit_preval(:,id), true );
[ ~, iy ] = min( xx.^2 + ( yy - 1 ).^2 );
fprintf( 'auc: %g index: %d\n', auc_, id )
fprintf( 'thr: %g, tpr: %g tnr: %g index: %d\n', tt(iy), yy(iy), 1 - xx(iy), iy )

nz = oo.glmnet_fit.beta(:,id) ~= 0;
yp = oo.fit_preval(:,id);

end

function [ Xs X ix ] = batch_filter( data, batch, repl )
%
% Batch correction using Bayesian fixed effects model and robust linear
% regression, filter probes using variance and entropy filter.
%

% batch effects correction
e = mdl_batch( double( data ), batch );
[ xpr0 xpr1 ] = deal( e(:,batch == 0), e(:,batch == 1) );
[ rep0 rep1 ] = deal( repl(batch == 0), repl(batch == 1) );

[ ~, ia ib ] = intersect( rep0, rep1 ); 

warning( 'off', 'stats:statrobustfit:IterationLimit' )
b = zeros( size( xpr0, 1 ), 2 );
for i=1:size( xpr0, 1 )
	[ b(i,:) ] = robustfit( xpr0(i,ia), xpr1(i,ib), 'huber' );
end
warning( 'on', 'stats:statrobustfit:IterationLimit' )
xpr0 = bsxfun( @plus, bsxfun( @times, xpr0, b(:,2) ), b(:,1) );

% drop replicates
xpr0(:,ia) = [];

X = [ xpr0 xpr1 ];

% filter probes
thr = 60;
i1 = geneentropyfilter( X, 'Percentile', thr );
i2 = genevarfilter( X, 'Percentile', thr );
Xs = X(i1 & i2,:);
ix = i1 & i2;

end
