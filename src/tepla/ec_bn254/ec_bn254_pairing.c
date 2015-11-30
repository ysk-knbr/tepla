//=============================================================================
// pairing over BN curve (ec_bn254_pairing) implementation with GMP
//=============================================================================
//  2015.10.31 created by kanbara
//============================================================================

#include "ec_bn254_lcl.h"

#define xcoord(p)   (p->x)
#define ycoord(p)   (p->y)
#define zcoord(p)   (p->z)

#define field(p)   (p->ec->field)
#define curve(p)   (p->ec)

//-------------------------------------------
//  precomputation for pairing
//-------------------------------------------
void ec_bn254_pairing_precomp_beuchat(EC_PAIRING p)
{
	pairing_precomp_p precomp;

	int s[] = {0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,-1,-1,0,0,0,0,0,0,1,1};
	int t[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,1};

	int *sbuff, *tbuff;

	precomp = (pairing_precomp_p)malloc(sizeof(struct ec_pairing_precomp_st));

	sbuff = (int *)malloc(sizeof(s));
	memcpy(sbuff, s, sizeof(s));

	precomp->si = sbuff;
	precomp->slen = sizeof(s)/sizeof(*s);

	tbuff = (int *)malloc(sizeof(t));
	memcpy(tbuff, t, sizeof(t));

	precomp->ti = tbuff;
	precomp->tlen = sizeof(t)/sizeof(*t);

	p->precomp = (void*)precomp;
}

void ec_bn254_pairing_precomp_aranha(EC_PAIRING p)
{
	pairing_precomp_p precomp;

	int s[] = {0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,1,1}; //65
	int t[] = {-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,-1}; //63

	int *sbuff, *tbuff;

	precomp = (pairing_precomp_p)malloc(sizeof(struct ec_pairing_precomp_st));

	sbuff = (int *)malloc(sizeof(s));
	memcpy(sbuff, s, sizeof(s));

	precomp->si = sbuff;
	precomp->slen = sizeof(s)/sizeof(*s);

	tbuff = (int *)malloc(sizeof(t));
	memcpy(tbuff, t, sizeof(t));

	precomp->ti = tbuff;
	precomp->tlen = sizeof(t)/sizeof(*t);

	p->precomp = (void*)precomp;
}

//-------------------------------------------
//  pairing (Beuchat)
//  point doubling
//-------------------------------------------
void ec_bn254_pairing_dob_beuchat(EC_POINT T, Element l0, Element l2, Element l4, const EC_POINT P)
{
	Element *t = field(T)->tmp;

	bn254_fp2_sqr(t[4], zcoord(T)); //ZT_2 = ZT^2
	bn254_fp2_sqr(l0, xcoord(T));   //l0 = XT^2

	bn254_fp2_sqr(t[0], ycoord(T)); //tmp0 = YT^0
	bn254_fp2_sqr(t[1], t[0]);      //tmp1 = tmp0^2

	bn254_fp2_add(l2, xcoord(T), t[0]); //l2 = XT + tmp0
	bn254_fp2_sqr(l2, l2);       //l2 = l2^3
	bn254_fp2_sub(l2, l2, l0);   //l2 = l2 - l0
	bn254_fp2_sub(l2, l2, t[1]); //l2 = l2 - tmp1
	bn254_fp2_add(l2, l2, l2);   //l2 = 2*l2

	bn254_fp2_add(t[2], l0, l0);   //tmp2 = 3*l0
	bn254_fp2_add(t[2], t[2], l0); //tmp2 = 3*l0
	bn254_fp2_sqr(t[3], t[2]);     //tmp3 = tmp2^2

	bn254_fp2_add(l4, t[2], xcoord(T)); //l4 = tmp2 + XT
	bn254_fp2_sqr(l4, l4);       //l4 = l4^2
	bn254_fp2_sub(l4, l4, l0);   //l4 = l4 - l0
	bn254_fp2_sub(l4, l4, t[3]); //l4 = l4 - tmp3
	bn254_fp2_sub(xcoord(T), t[3], l2); //XT = tmp3 - l2
	bn254_fp2_sub(xcoord(T), xcoord(T), l2); //XT = XT - l2

	bn254_fp2_add(zcoord(T), ycoord(T), zcoord(T)); //ZT = YT + ZT
	bn254_fp2_sqr(zcoord(T), zcoord(T));       //ZT = ZT^2
	bn254_fp2_sub(zcoord(T), zcoord(T), t[0]); //ZT = ZT - tmp0
	bn254_fp2_sub(zcoord(T), zcoord(T), t[4]); //ZT = ZT - ZT_2

	bn254_fp2_sub(ycoord(T), l2, xcoord(T));   //YT = l2 - XT
	bn254_fp2_mul(ycoord(T), ycoord(T), t[2]); //YT = YT*tmp2
	bn254_fp2_add(t[1], t[1], t[1]);   //tmp1 = 2*tmp1
	bn254_fp2_add(t[1], t[1], t[1]);   //tmp1 = 2*tmp1
	bn254_fp2_add(t[1], t[1], t[1]);   //tmp1 = 2*tmp1
	bn254_fp2_sub(ycoord(T), ycoord(T), t[1]); //YT = YT - tmp1

	bn254_fp2_mul(l2, t[2], t[4]);      //l2 = tmp2*ZT
	bn254_fp2_mul(l0, zcoord(T), t[4]); //l0 = ZT*ZT_2

	bn254_fp2_neg(l2, l2);   //l2 = -l2
	bn254_fp2_mul_p(l2, l2, xcoord(P)); //l2 = l2*XP

	bn254_fp2_add(t[1], t[0], t[0]); //tmp1 = 2*tmp0
	bn254_fp2_add(t[1], t[1], t[1]); //tmp1 = 2*tmp1
	bn254_fp2_sub(l4, l4, t[1]);     //l4 = l4 - tmp1

	bn254_fp2_mul_p(l0, l0, ycoord(P)); //l0 = l0*YP

	bn254_fp2_add(l0, l0, l0);
	bn254_fp2_add(l2, l2, l2);
}

//-------------------------------------------
//  pairing (Aranha)
//  point doubling in Jacobian coordinates
//-------------------------------------------
void ec_bn254_pairing_dob_aranha(EC_POINT T, Element l0, Element l2, Element l4, const EC_POINT P)
{
	//printf("dob\n");
	Element *t = field(T)->tmp;
	Element *k = field(P)->tmp;
	// 2T = R
 	// calculate R and l
 	bn254_fp2_sqr(t[0], xcoord(T)); 			// t0 = Tx^2
 	bn254_fp2_sqr(t[2], zcoord(T)); 			// t2 = Tz^2
 	bn254_fp2_dob(t[1], t[0]);      			// t1 = 2*t0 
 	bn254_fp2_mul(zcoord(T), ycoord(T), zcoord(T)); // Tz = Ty*Tz
 	bn254_fp2_add(t[0], t[0], t[1]); 			// t0 = t0+t1
 	bn254_fp2_sqr(t[3], ycoord(T));  			// t3 = Ty^2
 	bn254_fp2_div_2(t[0], t[0]);				// t0 = t0/2
 	bn254_fp2_mul(t[1], t[0], t[2]); 			// t1 = t0*t2
 	bn254_fp2_mul(t[4], t[0], xcoord(T)); 		// t4 = t0*Tx 
 	
 	bn254_fp_neg(k[0], xcoord(P)); 				// k0 = -Px
 	bn254_fp2_mul_p(l2, t[1], k[0]);	 	// l2 = l_(1,0) = t1*k0
 	//bn254_fp2_mul_p(l2, t[1], xcoord(P));	 	// l2 = l_(1,0) = t1*k0
 	//bn254_fp2_neg(l2,l2);
 	
 	bn254_fp2_sub(l4, t[4], t[3]);        		// l4 = l_(1,1) = t4-t3
 	bn254_fp2_mul(t[2], zcoord(T), t[2]); 		// t2 = Rz*t2
 	bn254_fp2_mul(t[1], t[3], xcoord(T)); 		// t1 = t3*Tx
 	bn254_fp2_mul_p(l0, t[2], ycoord(P));   	// l0 = l_(0,0) = t2*Py
 	bn254_fp2_dob(ycoord(T), t[1]); 	  		// Ty = 2*t1
 	bn254_fp2_sqr(xcoord(T), t[0]);       		// Tx = t0^2
 	bn254_fp2_sub(xcoord(T), xcoord(T), ycoord(T)); // Rx = Tx-Ty
 	bn254_fp2_sub(t[1], t[1], xcoord(T)); 		// t1 = t1-Rx
 	
 	bn254_fp2_muln(t[2], t[3], t[3]); 			// t2 = t3^2
 	bn254_fp2_OP1_2(t[2], t[2]);
 	bn254_fp2_muln(t[1], t[0], t[1]);			// t1 = t0*t1
 	bn254_fp2_OP1_2(t[1], t[1]);
 	bn254_fp2_subn(t[1], t[1], t[2]); 			// t1 = t1-t2
 	bn254_fp2_OP2(t[1], t[1]);
 	bn254_fp2_mod(ycoord(T), t[1]); 			// Ry = t1 mod p
}

//-------------------------------------------
//  pairing (Beuchat)
//  point addition 
//-------------------------------------------
void ec_bn254_pairing_add_beuchat(EC_POINT T, Element l0, Element l2, Element l4, const EC_POINT Q, const EC_POINT P)
{
	Element *t = field(T)->tmp;

	bn254_fp2_sqr(t[5], zcoord(T));   //ZT2 = ZT^2
	bn254_fp2_sqr(t[6], ycoord(Q));   //YQ2 = YQ^2

	bn254_fp2_mul(t[0], xcoord(Q), t[5]);    //t0 = XQ*ZT^2
	bn254_fp2_add(l2, ycoord(Q), zcoord(T)); //l2 = YQ*ZT
	bn254_fp2_sqr(l2, l2);         //l2 = l2^2
	bn254_fp2_sub(l2, l2, t[6]);   //l2 = l2 - YQ^2
	bn254_fp2_sub(l2, l2, t[5]);   //l2 = l2 - ZT^2
	bn254_fp2_mul(l2, l2, t[5]);   //l2 = l2*ZT^2

	bn254_fp2_sub(t[0], t[0], xcoord(T)); //tmp0 = tmp0 - XT
	bn254_fp2_sqr(t[1], t[0]);            //tmp1 = tmp0^2

	bn254_fp2_add(t[2], t[1], t[1]); //tmp2 = 2*tmp1
	bn254_fp2_add(t[2], t[2], t[2]); //tmp2 = 2*tmp2

	bn254_fp2_mul(t[3], t[0], t[2]); //tmp3 = tmp0*tmp2

	bn254_fp2_sub(t[4], l2, ycoord(T));   //tmp4 = l2 - YT
	bn254_fp2_sub(t[4], ycoord(T), t[4]); //tmp4 = YT - tmp4

	bn254_fp2_mul(l4, t[4], xcoord(Q));   //l4 = tmp4*XQ
	bn254_fp2_mul(t[2], t[2], xcoord(T)); //tmp2 = tmp2*XT
	bn254_fp2_sqr(xcoord(T), t[4]);       //XT = tmp4^2
	bn254_fp2_sub(xcoord(T), xcoord(T), t[3]); //XT = XT - t3
	bn254_fp2_sub(xcoord(T), xcoord(T), t[2]); //XT = XT - t2
	bn254_fp2_sub(xcoord(T), xcoord(T), t[2]); //XT = XT - t2

	bn254_fp2_add(zcoord(T), zcoord(T), t[0]); //ZT = ZT + tmp0
	bn254_fp2_sqr(zcoord(T), zcoord(T));       //ZT = ZT^2
	bn254_fp2_sub(zcoord(T), zcoord(T), t[5]); //ZT = ZT - ZT^2
	bn254_fp2_sub(zcoord(T), zcoord(T), t[1]); //ZT = ZT - tmp1

	bn254_fp2_sqr(t[5], zcoord(T));
	bn254_fp2_add(t[1], ycoord(Q), zcoord(T)); //tmp1 = YQ + ZT
	bn254_fp2_sub(t[2], xcoord(T), t[2]);      //tmp2 = XT - tmp2
	bn254_fp2_mul(t[2], t[2], t[4]);           //tmp2 = tmp2*tmp4

	bn254_fp2_mul(t[0], ycoord(T), t[3]); //tmp0 = YT*tmp3
	bn254_fp2_add(t[0], t[0], t[0]);      //tmp0 = 2*tmp0
	bn254_fp2_sub(ycoord(T), t[2], t[0]); //YT = tmp2 - tmp0

	bn254_fp2_sqr(t[1], t[1]);       //tmp1 = tmp1^2
	bn254_fp2_sub(t[1], t[1], t[6]); //tmp1 = tmp1 - YQ^2
	bn254_fp2_sub(t[1], t[5], t[1]); //tmp1 = ZT^2 - tmp1

	bn254_fp2_add(l4, l4, l4);   //l4 = 2*l4
	bn254_fp2_sub(l4, t[1], l4); //l4 = tmp1 - l4

	bn254_fp2_mul_p(l0, zcoord(T), ycoord(P)); // l0 = ZT * YP
	bn254_fp2_mul_p(l2, t[4], xcoord(P));      // l2 = tmp4 * XP

	bn254_fp2_add(l0, l0, l0);
	bn254_fp2_add(l2, l2, l2);
}

//-------------------------------------------
//  pairing (Aranha)
//  point addition in Jacobian coordinates
//-------------------------------------------
void ec_bn254_pairing_add_aranha(EC_POINT T, Element l0, Element l2, Element l4, const EC_POINT Q, const EC_POINT P)
{
	Element *t = field(T)->tmp;
	Element *k = field(P)->tmp;
	//printf("add\n");

	// T + Q = R
	bn254_fp2_sqr(t[1], zcoord(T));        		// t1 = Tz^2
	bn254_fp2_mul(t[3], xcoord(Q), t[1]); 		// t3 = Qx*t1
	bn254_fp2_mul(t[1], t[1], zcoord(T));  		// t1 = t1*Tz
	bn254_fp2_sub(t[3], t[3], xcoord(T));  		// t3 = t3-Tx
	bn254_fp2_mul(t[4], t[1], ycoord(Q));  		// t4 = t1*Qy
	bn254_fp2_mul(zcoord(T), zcoord(T), t[3]); 	// Rz = Z1*t3
	bn254_fp2_sub(t[0], t[4], ycoord(T));  		// t0 = t4-Ty
	bn254_fp2_sqr(t[1], t[3]);             		// t1 = t3^2
	bn254_fp2_mul(t[4], t[1], t[3]);	   		// t4 = t1*t3
	bn254_fp2_mul(t[1], t[1], xcoord(T));  		// t1 = t1*Tx
	bn254_fp2_sqr(xcoord(T), t[0]);		   		// Rx = t0^2
	bn254_fp2_dob(t[3], t[1]);			   		// t3 = 2*t1
	bn254_fp2_sub(xcoord(T), xcoord(T), t[3]); 	// Rx = Rx-t3
	bn254_fp2_sub(xcoord(T), xcoord(T), t[4]); 	// Rx = Rx-t4
	bn254_fp2_sub(t[1], t[1], xcoord(T));  		// t1 = t1-Rx
	bn254_fp2_muln(t[2], t[0], t[1]); 			// t2 = t0*t1
	bn254_fp2_OP1_2(t[2], t[2]);
	bn254_fp2_muln(t[3], t[4], ycoord(T)); 		// t3 = t4*Ty
	bn254_fp2_OP1_2(t[3], t[3]);
	bn254_fp2_sub(t[2], t[2], t[3]); 			// t2 = t2-t3
	bn254_fp2_OP2(t[2], t[2]);
	bn254_fp2_mod(ycoord(T), t[2]);  			// Ry = t2 mod p
	
	bn254_fp_neg(k[0], xcoord(P));
	bn254_fp2_mul_p(l2, t[0], k[0]); 		// l2 = l_(1,0) = t0*Px
	//bn254_fp2_mul_p(l2, t[0], xcoord(P)); 		// l2 = l_(1,0) = t0*Px	
	//bn254_fp2_neg(l2,l2);
	
	bn254_fp2_muln(t[2], t[0], xcoord(Q)); 		// t2 = t0*Qx
	bn254_fp2_OP1_2(t[2], t[2]);
	bn254_fp2_muln(t[3], zcoord(T), ycoord(Q)); // t3 = Rz*Qy
	bn254_fp2_OP1_2(t[3], t[3]);
	bn254_fp2_sub(t[2], t[2], t[3]); 			// t2 = t2-t3
	bn254_fp2_OP2(t[2], t[2]);
	bn254_fp2_mod(l4, t[2]);					// l4 = l_(1,1) = t2 mod p
	bn254_fp2_mul_p(l0, zcoord(T), ycoord(P)); 	// l0 = l_(0,0) = Rz*Py
}

//-------------------------------------------
//  pairing (Beuchat)
//  miller's algorithm
//-------------------------------------------
void ec_bn254_pairing_miller_beuchat(Element z, const EC_POINT Q, const EC_POINT P, const EC_PAIRING p)
{
	int len, *s, i;

	EC_POINT T, R, S;

	Element xq, yq;
	Element f, l0, l2, l4;

	//--------------------------------
	//   init
	//--------------------------------
	element_init(f, z->field);
	element_init(xq, field(Q));
	element_init(yq, field(Q));

	point_init(R, curve(Q));
	point_init(T, curve(Q));
	point_init(S, curve(Q));

	element_init(l0, field(Q));
	element_init(l2, field(Q));
	element_init(l4, field(Q));

	len = ((pairing_precomp_p)(p->precomp))->slen;  // s = PAIRING->precomp->si
	s = ((pairing_precomp_p)(p->precomp))->si;

	bn254_fp12_set_one(f);        // f = 1
	ec_bn254_fp2_point_set(T, Q); // T = Q
	ec_bn254_fp2_neg(R, Q);       // R = -Q

	//-------------------------------
	//  Miller loop
	//-------------------------------
	for (i =len-2; i >= 0; i--)
	{
		ec_bn254_pairing_dob_beuchat(T, l0, l2, l4, P);   //T = 2T, l = l(P)

		bn254_fp12_sqr(f, f);             // f = f^2*l
		bn254_fp12_mul_L(f, l0, l2, l4);  //

		if ( s[i] )
		{
			if ( s[i] < 0 )
			{
				ec_bn254_pairing_add_beuchat(T, l0, l2, l4, R, P);   // T = T - Q, l = l(P)
			}
			else
			{
				ec_bn254_pairing_add_beuchat(T, l0, l2, l4, Q, P);   // T = T + Q, l = l(P)
			}

			bn254_fp12_mul_L(f, l0, l2, l4);
		}
	}

	//--------------------------------
	//   addition part
	//--------------------------------
	ec_bn254_tw_frob(S, Q);
	ec_bn254_pairing_add_beuchat(T, l0, l2, l4, S, P);   //addtion part 2
	bn254_fp12_mul_L(f, l0, l2, l4);

	ec_bn254_tw_frob2(S, Q);
	ec_bn254_fp2_neg(S, S);
	ec_bn254_pairing_add_beuchat(T, l0, l2, l4, S, P);   //addtion part 2
	bn254_fp12_mul_L(f, l0, l2, l4);

	bn254_fp12_set(z, f);

	//--------------------------------
	//   relase
	//--------------------------------
	element_clear(f);
	element_clear(l0);
	element_clear(l2);
	element_clear(l4);
	element_clear(xq);
	element_clear(yq);

	point_clear(T);
	point_clear(R);
	point_clear(S);
}

//-------------------------------------------
//  pairing (Aranha)
//  miller's algorithm
//-------------------------------------------
void ec_bn254_pairing_miller_aranha(Element z, const EC_POINT Q, const EC_POINT P, const EC_PAIRING p)
{
	Element d, e, f, l0, l2, l4;

	int len, *s, i;

	EC_POINT T, S;

	//--------------------------------
	//   init
	//--------------------------------
	element_init(d, z->field);
	element_init(e, z->field);
	element_init(f, z->field);
	
	point_init(T, curve(Q));
	point_init(S, curve(Q));

	element_init(l0, field(Q));
	element_init(l2, field(Q));
	element_init(l4, field(Q));

	len = ((pairing_precomp_p)(p->precomp))->slen-1;  // s = PAIRING->precomp->si
	s = ((pairing_precomp_p)(p->precomp))->si;

	ec_bn254_fp2_point_set(T, Q);
	bn254_fp12_set_one(d);
	bn254_fp12_set_one(e);

	ec_bn254_pairing_dob_aranha(T, l0, l2, l4, P); // T = 2Q, l = l(P) 		
	bn254_fp12_mul_L(d, l0, l2, l4);  

	if( s[len-1] )
	{ 
		ec_bn254_pairing_add_aranha(T, l0, l2, l4, Q, P); // T = T+Q
		bn254_fp12_mul_L(e, l0, l2, l4);
	}

	bn254_fp12_mul(f, d, e);

	for(i = len-2 ; i >= 0 ; i--)
	{
		ec_bn254_pairing_dob_aranha(T, l0, l2, l4, P);   // T = 2T
		bn254_fp12_sqr(f, f);             		  // f = f^2*l
		bn254_fp12_mul_L(f, l0, l2, l4);  

		if( s[i] )
		{
			ec_bn254_pairing_add_aranha(T, l0, l2, l4, Q, P); // T = T+Q
			bn254_fp12_mul_L(f, l0, l2, l4);  
		}
	}

	ec_bn254_fp2_neg(T, T); // T = -T
   	bn254_fp12_conj(f, f);  // f = f^(p^6)
   
   	bn254_fp12_set_one(d);
	bn254_fp12_set_one(e);

    ec_bn254_tw_frob(S, Q);
	ec_bn254_pairing_add_aranha(T, l0, l2, l4, S, P);   //addtion part 
	bn254_fp12_mul_L(d, l0, l2, l4);

	ec_bn254_tw_frob2(S, Q);
	ec_bn254_fp2_neg(S, S);
	ec_bn254_pairing_add_aranha(T, l0, l2, l4, S, P);   //addtion part 
	bn254_fp12_mul_L(e, l0, l2, l4);

	bn254_fp12_mul(d, d, e); // d = d*e
	bn254_fp12_mul(z, f, d); // f = f*d = f*(d*e)
	
	//--------------------------------
	//   relase
	//--------------------------------
	element_clear(d);
	element_clear(e);
	element_clear(f);
	element_clear(l0);
	element_clear(l2);
	element_clear(l4);
	point_clear(T);
	point_clear(S);
}

void ec_bn254_pairing_finalexp(Element z, const Element x, const EC_PAIRING p)
{
	Element *t = z->field->tmp;
	int len, *u;
	len = ((pairing_precomp_p)(p->precomp))->tlen;
	u = ((pairing_precomp_p)(p->precomp))->ti;

	//------------------------------------------------------------
	// (p^4-p^2+1)/r = lambda3*p^3+lambda2*p^2+lambda1*p+lambda0
	// t := -(2^62+2^55+1)
	// lambda3 := 1
	// lambda2 := 6t^2+1
	// lambda1 := -36t^3-18t^2-12t+1
	// lambda0 := -36t^3-30t^2-18t-2
	//------------------------------------------------------------

	//------------------------------------------------------------
	//	calculate x^(p^6-1)
	//------------------------------------------------------------
	bn254_fp12_conj(t[0], x);     	// t0  = conjugate of x
	bn254_fp12_inv(t[1], x);      	// t1  = x^-1
	bn254_fp12_mul(z, t[0], t[1]); 	// z  = t0*t1
								   	// z  = x^(p^6-1)
	//------------------------------------------------------------
	//	calculate x^(p^6-1)(p^2+1)
	//------------------------------------------------------------
	bn254_fp12_frob_p2(t[0], z);   	// t0 = z^(p^2)
	bn254_fp12_mul(z, z, t[0]);    	// z  = z*t0
								   	// z  = x^{(p^6-1)(p^2+1)}
	//------------------------------------------------------------
	//	calculate x^{(p^6-1)(p^2+1)(p^4-p^2+1)/r}
	//------------------------------------------------------------
	if(p->type == Pairing_ECBN254a)
	{
		bn254_fp12_pow_forpairing_beuchat(t[7], z, u, len);		// t7 = z^t
		bn254_fp12_pow_forpairing_beuchat(t[8], t[7], u, len);	// t8 = z^(t^2)
		bn254_fp12_pow_forpairing_beuchat(t[9], t[8], u, len);	// t9 = z^(t^3)
	}	

	if(p->type == Pairing_ECBN254b)
	{
		bn254_fp12_pow_forpairing_karabina(t[7], z, u, len);	// t7 = z^t
		bn254_fp12_pow_forpairing_karabina(t[8], t[7], u, len);	// t8 = z^(t^2)
		bn254_fp12_pow_forpairing_karabina(t[9], t[8], u, len);	// t9 = z^(t^3)
	}

	bn254_fp12_frob_p(t[0], z);			// t0 = z^p
	bn254_fp12_frob_p2(t[1], z); 		// t1 = z^(p^2)
	bn254_fp12_frob_p3(t[2], z); 		// t2 = z^(p^3)

	bn254_fp12_mul(t[0], t[0], t[1]); 	// t0 = z^p * z^(p^2)
	bn254_fp12_mul(t[0], t[0], t[2]); 	// t0 = z^p * z^(p^2) * z^(p^3)

	bn254_fp12_conj(t[1], z);			// t1 = 1/z

	bn254_fp12_frob_p2(t[2], t[8]); 	// t2 = (z^(t^2))^p^2

	bn254_fp12_frob_p(t[3], t[7]); 		// t3 = (z^x)^p
	bn254_fp12_conj(t[3], t[3]);		// t3 = 1/(z^x)^p

	bn254_fp12_frob_p(t[4], t[8]); 		// t4 = (z^(t^2))^p
	bn254_fp12_mul(t[4], t[4], t[7]); 	// t4 = z^t * (z^(t^2))^p
	bn254_fp12_conj(t[4], t[4]);		// t4 = 1/(z^t * (z^(t^2))^p)

	bn254_fp12_conj(t[5], t[8]);		// t5 = 1/(z^(t^2))

	bn254_fp12_frob_p(t[6], t[9]); 		// t6 = (z^(t^3))^p
	bn254_fp12_mul(t[6], t[6], t[9]); 	// t6 = z^(t^3) * (z^(t^3))^p
	bn254_fp12_conj(t[6], t[6]);		// t6 = 1/(z^(t^3) * (z^(t^3))^p)

	bn254_fp12_sqr(t[7], t[6]);			// t7 = t6^2
	bn254_fp12_mul(t[7], t[7], t[4]);	// t7 = t7 * t4
	bn254_fp12_mul(t[7], t[7], t[5]);	// t7 = t7 * t5
	bn254_fp12_mul(t[8], t[3], t[5]);	// t8 = t3 * t5
	bn254_fp12_mul(t[8], t[8], t[7]);	// t8 = t8 * t7
	bn254_fp12_mul(t[7], t[7], t[2]); 	// t7 = t7 * t2
	bn254_fp12_sqr(t[8], t[8]);			// t8 = t8^2
	bn254_fp12_mul(t[8], t[8], t[7]);	// t8 = t8 * t7
	bn254_fp12_sqr(t[8], t[8]);			// t8 = t8^2
	bn254_fp12_mul(t[7], t[8], t[1]);	// t7 = t8 * t1
	bn254_fp12_mul(t[8], t[8], t[0]);	// t8 = t8 * t0
	bn254_fp12_sqr(t[7], t[7]); 		// t7 = t7^2
	bn254_fp12_mul(z, t[7], t[8]);		// t7 = t7 * t8

}

void ec_bn254_pairing_aranha(Element z, const EC_POINT Q, const EC_POINT P, const EC_PAIRING p)
{
	ec_bn254_pairing_miller_aranha(z, Q, P, p);
	ec_bn254_pairing_finalexp(z, z, p);
}

void ec_bn254_pairing_beuchat(Element z, const EC_POINT Q, const EC_POINT P, const EC_PAIRING p)
{
	ec_bn254_pairing_miller_beuchat(z, Q, P, p);
	ec_bn254_pairing_finalexp(z, z, p);
}

void ec_bn254_double_pairing_beuchat(Element z, const EC_POINT Q1, const EC_POINT P1, const EC_POINT Q2, const EC_POINT P2, const EC_PAIRING p)
{
	Element z1, z2;
	element_init(z1, p->g3);
	element_init(z2, p->g3);

	ec_bn254_pairing_miller_beuchat(z1, Q1, P1, p);
	ec_bn254_pairing_miller_beuchat(z2, Q2, P2, p);
	bn254_fp12_mul(z, z1, z2);
	ec_bn254_pairing_finalexp(z, z, p);

	element_clear(z1);
	element_clear(z2);
}


void ec_bn254_double_pairing_aranha(Element z, const EC_POINT Q1, const EC_POINT P1, const EC_POINT Q2, const EC_POINT P2, const EC_PAIRING p)
{
	Element z1, z2;
	element_init(z1, p->g3);
	element_init(z2, p->g3);

	ec_bn254_pairing_miller_aranha(z1, Q1, P1, p);
	ec_bn254_pairing_miller_aranha(z2, Q2, P2, p);
	bn254_fp12_mul(z, z1, z2);
	ec_bn254_pairing_finalexp(z, z, p);

	element_clear(z1);
	element_clear(z2);
}
