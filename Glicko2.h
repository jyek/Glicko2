/*******************************************************************************
 GLICKO 2
 Implements the Glicko 2 ranking system
 By Justin Yek (modified from work by Heungsub Lee)
 *******************************************************************************/

@class Glicko2;

@interface Glicko2 : NSObject

// Ratings
@property float mu;        // mean
@property float phi;       //
@property float sigma;     // uncertainty

// Other parameters
@property float alpha;
@property float diff;
@property float variance;

- (id) init;
- (id) initWithRatingMu:(float)_mu phi:(float)_phi sigma:(float)_sigma;
- (void) printRating;
- (void) scaleDown;
- (void) scaleUp;
- (float) reduceImpact;
- (float) determineSigma;
- (float) f:(float)x;
- (void) rateSingle:(float)score ratings:(Glicko2*)otherRating;
- (void) rate:(NSArray*)scores ratings:(NSArray*)ratings;

@end

