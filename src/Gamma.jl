# Table[Gamma[(2n+1)/2], {n,-170,171}]
const GAMMA_HALF = (
    5.6482208842233255e-306, -9.5737343987585367e-304,
    1.6131742461908134e-301, -2.7020668623696125e-299,  4.4989413258454048e-297,
   -7.4457478942741450e-295,  1.2248255286080968e-292, -2.0025897392742383e-290,
    3.2542083263206373e-288, -5.2555464470078292e-286,  8.4351520474475659e-284,
   -1.3454067515678868e-281,  2.1324697012351005e-279, -3.3586397794452833e-277,
    5.2562712548318684e-275, -8.1735018012635554e-273,  1.2628060282952193e-270,
   -1.9384072534331616e-268,  2.9560710614855715e-266, -4.4784476581506408e-264,
    6.7400637255167144e-262, -1.0076395269647488e-259,  1.4963446975426520e-257,
   -2.2071084288754117e-255,  3.2334138483024781e-253, -4.7046171492801056e-251,
    6.7981717807097526e-249, -9.7553765053184950e-247,  1.3901411520078855e-244,
   -1.9670497300911580e-242,  2.7637048707780770e-240, -3.8553682947354175e-238,
    5.3396850882085532e-236, -7.3420669962867607e-234,  1.0021921449931428e-231,
   -1.3579703564657085e-229,  1.8264701294463780e-227, -2.4383376228109146e-225,
    3.2307973502244618e-223, -4.2484985155451673e-221,  5.5442905627864434e-219,
   -7.1798562788084442e-217,  9.2261153182688507e-215, -1.1763297030792785e-212,
    1.4880570743952873e-210, -1.8675116283660855e-208,  2.3250519773157765e-206,
   -2.8714391919849839e-204,  3.5175130101816053e-202, -4.2737783073706505e-200,
    5.1499028603816338e-198, -6.1541339181560524e-196,  7.2926486930149221e-194,
   -8.5688622142925335e-192,  9.9827244796508015e-190, -1.1530046773996676e-187,
    1.3201903556226194e-185, -1.4984160536316730e-183,  1.6857180603356321e-181,
   -1.8795756372742298e-179,  2.0769310791880239e-177, -2.2742395317108862e-175,
    2.4675498919063115e-173, -2.6526161337992849e-171,  2.8250361824962384e-169,
   -2.9804131725335315e-167,  3.1145317652975404e-165, -3.2235403770829544e-163,
    3.3041288865100282e-161, -3.3536908198076786e-159,  3.3704592739067170e-157,
   -3.3536069775371835e-155,  3.3033028728741257e-153, -3.2207203010522726e-151,
    3.1079950905154430e-149, -2.9681353114422481e-147,  2.8048878693129244e-145,
   -2.6225701578075843e-143,  2.4258773959720155e-141, -2.2196778173143942e-139,
    2.0088084246695268e-137, -1.7978835400792264e-135,  1.5911269329701154e-133,
   -1.3922360663488510e-131,  1.2042841973917561e-129, -1.0296629887699515e-127,
    8.7006522551060899e-126, -7.2650446330135850e-124,  5.9936618222362076e-122,
   -4.8848343851225092e-120,  3.9322916800236199e-118, -3.1261718856187778e-116,
    2.4540449302107406e-114, -1.9018848209133240e-112,  1.4549418879986928e-110,
   -1.0984811254390131e-108,  8.1836843845206475e-107, -6.0150080226226759e-105,
    4.3608808164014401e-103, -3.1180297837270296e-101,  2.1982109975275559e-099,
   -1.5277566432816513e-097,  1.0465133006479312e-095, -7.0639647793735354e-094,
    4.6975365782834011e-092, -3.0768864587756277e-090,  1.9845917659102799e-088,
   -1.2602157713530277e-086,  7.8763485709564232e-085, -4.8439543711382003e-083,
    2.9305923945386112e-081, -1.7437024747504736e-079,  1.0200659477290271e-077,
   -5.8653791994419057e-076,  3.3139392476846767e-074, -1.8392362824649956e-072,
    1.0023837739434226e-070, -5.3627531905973109e-069,  2.8154454250635882e-067,
   -1.4499543939077479e-065,  7.3222696892341270e-064, -3.6245234961708929e-062,
    1.7578938956428830e-060, -8.3499960043036945e-059,  3.8827481420012179e-057,
   -1.7666504046105542e-055,  7.8615943005169660e-054, -3.4197935207248802e-052,
    1.4534122463080741e-050, -6.0316608221785075e-049,  2.4428226329822955e-047,
   -9.6491494002800673e-046,  3.7149225191078259e-044, -1.3930959446654347e-042,
    5.0848001980288367e-041, -1.8051040703002370e-039,  6.2276090425358178e-038,
   -2.0862490292494990e-036,  6.7803093450608716e-035, -2.1357974436941746e-033,
    6.5141822032672324e-032, -1.9216837499638336e-030,  5.4767986873969256e-029,
   -1.5061196390341546e-027,  3.9912170434405096e-026, -1.0177603460773299e-024,
    2.4935128478894584e-023, -5.8597551925402271e-022,  1.3184449183215511e-020,
   -2.8346565743913349e-019,  5.8110459775022365e-018, -1.1331539656129361e-016,
    2.0963348363839318e-015, -3.6685859636718807e-014,  6.0531668400586031e-013,
   -9.3824086020908348e-012,  1.3604492473031710e-010, -1.8366064838592809e-009,
    2.2957581048241011e-008, -2.6401218205477163e-007,  2.7721279115751021e-006,
   -2.6335215159963470e-005,  2.2384932885968950e-004, -1.6788699664476712e-003,
    1.0912654781909863e-002, -6.0019601300504246e-002,  2.7008820585226911e-001,
   -9.4530872048294188e-001,  2.3632718012073547e+000, -3.5449077018110321e+000,
    1.7724538509055160e+000,  8.8622692545275801e-001,  1.3293403881791370e+000,
    3.3233509704478426e+000,  1.1631728396567449e+001,  5.2342777784553520e+001,
    2.8788527781504436e+002,  1.8712543057977883e+003,  1.4034407293483413e+004,
    1.1929246199460901e+005,  1.1332783889487856e+006,  1.1899423083962248e+007,
    1.3684336546556586e+008,  1.7105420683195732e+009,  2.3092317922314238e+010,
    3.3483860987355646e+011,  5.1899984530401251e+012,  8.5634974475162064e+013,
    1.4986120533153361e+015,  2.7724322986333718e+016,  5.4062429823350750e+017,
    1.1082798113786904e+019,  2.3828015944641843e+020,  5.3613035875444147e+021,
    1.2599063430729375e+023,  3.0867705405286968e+024,  7.8712648783481768e+025,
    2.0858851927622669e+027,  5.7361842800962338e+028,  1.6348125198274266e+030,
    4.8226969334909086e+031,  1.4709225647147271e+033,  4.6334060788513904e+034,
    1.5058569756267019e+036,  5.0446208683494513e+037,  1.7403941995805607e+039,
    6.1783994085109905e+040,  2.2551157841065115e+042,  8.4566841903994183e+043,
    3.2558234133037760e+045,  1.2860502482549915e+047,  5.2085035054327157e+048,
    2.1615289547545770e+050,  9.1864980577069524e+051,  3.9961266551025243e+053,
    1.7782763615206233e+055,  8.0911574449188360e+056,  3.7623882118872587e+058,
    1.7871344006464479e+060,  8.6676018431352723e+061,  4.2904629123519598e+063,
    2.1666837707377397e+065,  1.1158421419299359e+067,  5.8581712451321637e+068,
    3.1341216161457076e+070,  1.7080962807994106e+072,  9.4799343584367290e+073,
    5.3561629125167519e+075,  3.0797936746971323e+077,  1.8016792996978224e+079,
    1.0719991833202043e+081,  6.4855950590872363e+082,  3.9886409613386503e+084,
    2.4929006008366564e+086,  1.5829918815312768e+088,  1.0210297635876736e+090,
    6.6877449514992618e+091,  4.4473503927470091e+093,  3.0019615151042312e+095,
    2.0563436378463983e+097,  1.4291588283032468e+099,  1.0075569739537890e+101,
    7.2040323637695915e+102,  5.2229234637329539e+104,  3.8388487458437211e+106,
    2.8599423156535722e+108,  2.1592564483184470e+110,  1.6518311829636120e+112,
    1.2801691667967993e+114,  1.0049327959354874e+116,  7.9892157276871251e+117,
    6.4313186607881357e+119,  5.2415247085423306e+121,  4.3242578845474227e+123,
    3.6107553335970980e+125,  3.0510882568895478e+127,  2.6086804596405634e+129,
    2.2565085975890873e+131,  1.9744450228904514e+133,  1.7473838452580495e+135,
    1.5639085415059543e+137,  1.4153372300628886e+139,  1.2950335655075431e+141,
    1.1979060480944774e+143,  1.1200421549683363e+145,  1.0584398364450778e+147,
    1.0108100438050493e+149,  9.7543169227187261e+150,  9.5104589996507580e+152,
    9.3678021146559966e+154,  9.3209631040827166e+156,  9.3675679196031302e+158,
    9.5080814383971771e+160,  9.7457834743571066e+162,  1.0086885895959605e+165,
    1.0540795761277788e+167,  1.1120539528148066e+169,  1.1843374597477690e+171,
    1.2731627692288517e+173,  1.3813816046133041e+175,  1.5126128570515680e+177,
    1.6714372070419826e+179,  1.8636524858518106e+181,  2.0966090465832869e+183,
    2.3796512678720307e+185,  2.7247007017134751e+187,  3.1470293104790638e+189,
    3.6662891467081093e+191,  4.3078897473820284e+193,  5.1048493506477037e+195,
    6.1002949740240059e+197,  7.3508554436989271e+199,  8.9312893640941964e+201,
    1.0940829471015391e+204,  1.3511924396704007e+206,  1.6822345873896489e+208,
    2.1112044071740094e+210,  2.6706735750751219e+212,  3.4051088082207804e+214,
    4.3755648185637028e+216,  5.6663564400399951e+218,  7.3945951542521937e+220,
    9.7238926278416347e+222,  1.2884157731890166e+225,  1.7200350572073372e+227,
    2.3134471519438685e+229,  3.1347208908839418e+231,  4.2788940160565805e+233,
    5.8834792720777982e+235,  8.1486187918277505e+237,  1.1367323214599712e+240,
    1.5971089116512595e+242,  2.2599091099865322e+244,  3.2203704817308084e+246,
    4.6212316412837101e+248,  6.6776797216549611e+250,  9.7160239950079684e+252,
    1.4233975152686674e+255,  2.0995113350212844e+257,  3.1177743325066073e+259,
    4.6610726270973779e+261,  7.0149143037815538e+263,  1.0627595170229054e+266,
    1.6207082634599307e+268,  2.4877871844109937e+270,  3.8436311999149852e+272,
    5.9768465158678020e+274,  9.3537647973331102e+276,  1.4732179555799648e+279,
    2.3350504595942443e+281,  3.7244054830528196e+283,  5.9776708002997755e+285,
    9.6539383424841375e+287,  1.5687649806536723e+290,  2.5649307433687543e+292,
    4.2193110728416008e+294,  6.9829598255528493e+296,  1.1626628109545494e+299,
    1.9474602083488703e+301,  3.2814704510678464e+303,  5.5620924145599996e+305,
    9.4833675668247993e+307
)

# Gamma function for half-integer arguments, gammahalf(n) = gamma(n/2)
function gammahalf(n::Integer)::Float64
    if iseven(n)
        if n < 0
            iseven(n) && throw(DomainError(n, "gammahalf not implemented even arguments n < 0"))
        else
            fac(n÷2)
        end
    else
        i = (n + 1)÷2 + 170
        if i < 1
            0.0
        elseif i <= length(GAMMA_HALF)
            GAMMA_HALF[i]
        else
            Inf
        end
    end
end