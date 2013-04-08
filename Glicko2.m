/*******************************************************************************
 GLICKO 2
 Implements the Glicko 2 ranking system
 By Justin Yek (modified from work by Heungsub Lee)
 *******************************************************************************/

#import "Glicko2.h"

#define win 1           // The actual score for win
#define draw 0.5        // The actual score for draw
#define loss 0          // The actual score for loss

const float defaultMu = 1500;      // default rating mu in Glicko-1
const float defaultPhi = 350;      // default rating deviation in Glicko-1
const float defaultSigma = 0.6;    // default rating sigma in Glicko-1
const float tau = 0.5;             // contrains the change in volatility over time, reasonable choices are between 0.3 and 1.2
const float epsilon = 0.000001;    // convergence tolerance
const float ratio = 173.7178;      // scale down and scale up factor in Glicko 2

@implementation Glicko2

@synthesize mu, phi, sigma, alpha, diff, variance;

- (id) init {
    return [self initWithRatingMu:defaultMu phi:defaultPhi sigma:defaultSigma];
}

- (id) initWithRatingMu:(float)_mu phi:(float)_phi sigma:(float)_sigma {
    if (self = [super init]){
        mu = _mu;
        phi = _phi;
        sigma = _sigma;
    }
    return self;
}

- (void) printRating {
    CCLOG(@"Mu = %f, Phi = %f, Sigma = %f",mu,phi,sigma);
}

- (void) scaleDown {
    mu = (mu - defaultMu) / ratio;
    phi = phi / ratio;
}

- (void) scaleUp {
    mu = mu * ratio + defaultMu;
    phi = phi * ratio;
}

- (float) reduceImpact {
    // The original form is `g(RD)`. This function reduces the impact of games as a function of an opponent's RD.
    return 1 / sqrt(1 + (3 * pow(phi,2)/pow(M_PI,2)));
}

// Returns new sigma
- (float) determineSigma {

    // 1. Let a = ln(s^2), and define f(x)
    alpha = log(pow(sigma,2));

    // 2. Set the initial values of the iterative algorithm.
    float a = alpha;
    float b;
    if (pow(diff,2) > (pow(phi,2) + variance) ) {
        b = log(pow(diff,2) - pow(phi,2) - variance);
    }
    else {
        float k = 1;
        b = alpha - k * sqrt(pow(tau,2));
        while ([self f:(alpha - k * sqrt(pow(tau,2)))] < 0) {
            k += 1;
            b = alpha - k * sqrt(pow(tau,2));
        }
    }
    
    // 3. Let fA = f(A) and f(B) = f(B)
    float f_a = [self f:a];
    float f_b = [self f:b];

    //NSLog(@"  Determining sigma: A = %f B = %f F_A = %f F_B = %f",a,b,f_a,f_b);
    
    // 4. While |B-A| > e, carry out the following steps.
    // (a) Let C = A + (A - B)fA / (fB-fA), and let fC = f(C).
    // (b) If fCfB < 0, then set A <- B and fA <- fB; otherwise, just set
    //     fA <- fA/2.
    // (c) Set B <- C and fB <- fC.
    // (d) Stop if |B-A| <= e. Repeat the above three steps otherwise.

    while (fabsf(b - a) > epsilon){
        float c = a + (a - b) * f_a / (f_b - f_a);
        float f_c = [self f:c];
        if (f_c * f_b < 0){
            a = b;
            f_a = f_b;
        }
        else {
            f_a = f_a / 2;
        }
        b = c;
        f_b = f_c;
        //NSLog(@"  Iteration: A = %f B = %f F_A = %f F_B = %f",a,b,f_a,f_b);
    }
    
    // 5. Once |B-A| <= e, set s' <- e^(A/2)
    return pow(exp(1),a/2);
}

// This function is twice the conditional log-posterior density of
// phi, and is the optimality criterion.
- (float) f:(float)x {
    float tmp = pow(phi,2) + variance + exp(x);
    float a = exp(x) * (pow(diff,2) - tmp) / (2 * pow(tmp,2));
    float b = (x - alpha) / pow(tau,2);
    return a - b;
}

// Rate a single score
- (void) rateSingle:(float)score ratings:(Glicko2*)otherRating {
    NSArray *scores = [NSArray arrayWithObject:[NSNumber numberWithFloat:score]];
    NSArray *ratings = [NSArray arrayWithObject:otherRating];
    [self rate:scores ratings:ratings];
}

// Rate multiple scores
- (void) rate:(NSArray*)scores ratings:(NSArray*)ratings {
    
    // Step 2. For each player, convert the rating and RD's onto the
    //         Glicko-2 scale.
    [self scaleDown];
    
    // Step 3. Compute the quantity v. This is the estimated variance of the
    //         team's/player's rating based only on game outcomes.
    // Step 4. Compute the quantity difference, the estimated improvement in
    //         rating by comparing the pre-period rating to the performance
    //         rating based only on game outcomes.
    diff = 0;
    variance = 0;
    float variance_inv = 0;
    int i = 0;
    
    for (Glicko2 *thisRating in ratings){
        float actualScore = [(NSNumber*)scores[i] intValue];
        [thisRating scaleDown];
        float imp = [thisRating reduceImpact];
        float expectedScore = 1 / (1 + exp(-imp * (mu - thisRating.mu)));
        variance_inv += pow(imp,2) * expectedScore * (1 - expectedScore);
        diff += imp * (actualScore - expectedScore);
        //NSLog(@"OTHER mu = %f, phi = %f, sigma = %f, g = %f, e = %f",thisRating.mu,thisRating.phi,thisRating.sigma,imp,expectedScore);
        i++;
    }
    
    variance = 1 / variance_inv;
    diff *= variance;
    
    // Step 5. Determine the new value, Sigma', ot the sigma. This
    //         computation requires iteration.
    sigma = [self determineSigma];
    
    // Step 6. Update the rating deviation to the new pre-rating period
    //         value, Phi*.
    float phi_star = sqrt(pow(phi,2) + pow(sigma,2));
    
    // Step 7. Update the rating and RD to the new values, Mu' and Phi'.
    phi = 1 / sqrt(1 / pow(phi_star,2) + 1 / variance);
    mu = mu + pow(phi,2) * (diff / variance);                       // note: diff / variance = imp * (actualScore - expectedScore)
    
    // Step 8. Convert ratings and RD's back to original scale.
    [self scaleUp];
}

@end
