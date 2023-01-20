/*******************************************************************************
 * Copyright (C) 2017-2023 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#include "IntegrationPlan.h"

arma::mat generate_points(const unsigned dimension, const std::pair<arma::vec, arma::vec>& points) {
    const auto& [location, weight] = points;

    if(1 == dimension) return join_rows(location, weight);

    arma::mat int_pts;

    const auto order = location.n_rows;

    auto counter = 0;

    if(2 == dimension) {
        int_pts.set_size(order * order, 3);
        for(auto i = 0llu; i < order; ++i)
            for(auto j = 0llu; j < order; ++j) {
                int_pts(counter, 0) = location(i);
                int_pts(counter, 1) = location(j);
                int_pts(counter++, 2) = weight(i) * weight(j);
            }
        return int_pts;
    }

    if(3 == dimension) {
        int_pts.set_size(order * order * order, 4);
        for(auto i = 0llu; i < order; ++i)
            for(auto j = 0llu; j < order; ++j)
                for(auto k = 0llu; k < order; ++k) {
                    int_pts(counter, 0) = location(i);
                    int_pts(counter, 1) = location(j);
                    int_pts(counter, 2) = location(k);
                    int_pts(counter++, 3) = weight(i) * weight(j) * weight(k);
                }
        return int_pts;
    }

    throw std::invalid_argument("not supported");
}

template<IntegrationType> std::pair<arma::vec, arma::vec> generate_seeds(unsigned) { throw std::invalid_argument("not supported"); }

template<> std::pair<arma::vec, arma::vec> generate_seeds<IntegrationType::GAUSS>(const unsigned order) {
    arma::vec PTL(order, arma::fill::none), PTW(order, arma::fill::none);

    switch(order) {
    case 1: {
        PTL(0) = 0.;

        PTW(0) = 2.;

        break;
    }
    case 2: {
        PTL(0) = -1. / sqrt(3.);

        PTW(0) = 1.;

        break;
    }
    case 3: {
        PTL(0) = -sqrt(.6);
        PTL(1) = 0.;

        PTW(0) = 5. / 9.;
        PTW(1) = 8. / 9.;

        break;
    }
    case 4: {
        constexpr auto TMPA = 3. / 7.;
        const auto TMPB = 2. / 7. * sqrt(1.2);
        const auto TMPC = sqrt(30.) / 36.;

        PTL(0) = -sqrt(TMPA + TMPB);
        PTL(1) = -sqrt(TMPA - TMPB);

        PTW(0) = -TMPC + .5;
        PTW(1) = TMPC + .5;

        break;
    }
    case 5: {
        const auto TMPA = 2. * sqrt(10. / 7.);
        const auto TMPB = 13. * sqrt(70.);

        PTL(0) = -sqrt(5. + TMPA) / 3.;
        PTL(1) = -sqrt(5. - TMPA) / 3.;
        PTL(2) = 0.;

        PTW(0) = (322. - TMPB) / 900.;
        PTW(1) = (322. + TMPB) / 900.;
        PTW(2) = 512. / 900.;

        break;
    }
    case 6: {
        PTL(0) = -.932469514203152028;
        PTL(1) = -.661209386466264514;
        PTL(2) = -.238619186083196909;

        PTW(0) = .171324492379170345;
        PTW(1) = .360761573048138608;
        PTW(2) = .46791393457269105;

        break;
    }
    case 7: {
        PTL(0) = -.949107912342758525;
        PTL(1) = -.74153118559939444;
        PTL(2) = -.405845151377397167;
        PTL(3) = .0;

        PTW(0) = .129484966168869693;
        PTW(1) = .279705391489276668;
        PTW(2) = .381830050505118945;
        PTW(3) = .417959183673469388;

        break;
    }
    case 8: {
        PTL(0) = -.960289856497536232;
        PTL(1) = -.79666647741362674;
        PTL(2) = -.525532409916328986;
        PTL(3) = -.183434642495649805;

        PTW(0) = .101228536290376259;
        PTW(1) = .22238103445337447;
        PTW(2) = .313706645877887287;
        PTW(3) = .362683783378361983;

        break;
    }
    case 9: {
        PTL(0) = -.96816023950762609;
        PTL(1) = -.836031107326635794;
        PTL(2) = -.613371432700590397;
        PTL(3) = -.32425342340380893;
        PTL(4) = .0;

        PTW(0) = .081274388361574412;
        PTW(1) = .1806481606948574;
        PTW(2) = .260610696402935462;
        PTW(3) = .31234707704000284;
        PTW(4) = .33023935500125976;

        break;
    }
    case 10: {
        PTL(0) = -.97390652851717172;
        PTL(1) = -.86506336668898451;
        PTL(2) = -.679409568299024406;
        PTL(3) = -.433395394129247191;
        PTL(4) = -.14887433898163121;

        PTW(0) = .066671344308688138;
        PTW(1) = .14945134915058059;
        PTW(2) = .219086362515982044;
        PTW(3) = .26926671930999636;
        PTW(4) = .29552422471475287;

        break;
    }
    case 11: {
        PTL(0) = -.978228658146056993;
        PTL(1) = -.887062599768095299;
        PTL(2) = -.730152005574049324;
        PTL(3) = -.519096129206811816;
        PTL(4) = -.26954315595234497;
        PTL(5) = .0;

        PTW(0) = .055668567116173666;
        PTW(1) = .125580369464904625;
        PTW(2) = .186290210927734251;
        PTW(3) = .23319376459199048;
        PTW(4) = .26280454451024666;
        PTW(5) = .27292508677790063;

        break;
    }
    case 12: {
        PTL(0) = -.981560634246719251;
        PTL(1) = -.904117256370474857;
        PTL(2) = -.769902674194304687;
        PTL(3) = -.587317954286617447;
        PTL(4) = -.367831498998180194;
        PTL(5) = -.125233408511468916;

        PTW(0) = .047175336386511827;
        PTW(1) = .10693932599531843;
        PTW(2) = .160078328543346226;
        PTW(3) = .203167426723065922;
        PTW(4) = .23349253653835481;
        PTW(5) = .249147045813402785;

        break;
    }
    case 13: {
        PTL(0) = -.98418305471858815;
        PTL(1) = -.917598399222977965;
        PTL(2) = -.801578090733309913;
        PTL(3) = -.642349339440340221;
        PTL(4) = -.448492751036446853;
        PTL(5) = -.230458315955134794;
        PTL(6) = 0.;

        PTW(0) = .04048400476531588;
        PTW(1) = .092121499837728448;
        PTW(2) = .138873510219787239;
        PTW(3) = .178145980761945738;
        PTW(4) = .207816047536888502;
        PTW(5) = .22628318026289724;
        PTW(6) = .23255155323087391;

        break;
    }
    case 14: {
        PTL(0) = -.986283808696812339;
        PTL(1) = -.928434883663573517;
        PTL(2) = -.827201315069764993;
        PTL(3) = -.68729290481168547;
        PTL(4) = -.515248636358154092;
        PTL(5) = -.31911236892788976;
        PTL(6) = -.108054948707343662;

        PTW(0) = .035119460331751863;
        PTW(1) = .08015808715976021;
        PTW(2) = .121518570687903185;
        PTW(3) = .157203167158193535;
        PTW(4) = .185538397477937814;
        PTW(5) = .205198463721295604;
        PTW(6) = .21526385346315779;

        break;
    }
    case 15: {
        PTL(0) = -.987992518020485429;
        PTL(1) = -.9372733924007059;
        PTL(2) = -.848206583410427216;
        PTL(3) = -.724417731360170047;
        PTL(4) = -.570972172608538848;
        PTL(5) = -.39415134707756337;
        PTL(6) = -.201194093997434522;
        PTL(7) = .0;

        PTW(0) = .0307532419961172684;
        PTW(1) = .0703660474881081247;
        PTW(2) = .107159220467171935;
        PTW(3) = .139570677926154314;
        PTW(4) = .166269205816993934;
        PTW(5) = .18616100001556221;
        PTW(6) = .198431485327111577;
        PTW(7) = .202578241925561273;

        break;
    }
    case 16: {
        PTL(0) = -.989400934991649933;
        PTL(1) = -.944575023073232576;
        PTL(2) = -.865631202387831744;
        PTL(3) = -.75540440835500303;
        PTL(4) = -.61787624440264375;
        PTL(5) = -.458016777657227386;
        PTL(6) = -.281603550779258913;
        PTL(7) = -.09501250983763744;

        PTW(0) = .027152459411754095;
        PTW(1) = .062253523938647893;
        PTW(2) = .095158511682492785;
        PTW(3) = .12462897125553387;
        PTW(4) = .149595988816576732;
        PTW(5) = .169156519395002538;
        PTW(6) = .18260341504492359;
        PTW(7) = .1894506104550685;

        break;
    }
    case 17: {
        PTL(0) = -.990575475314417336;
        PTL(1) = -.950675521768767761;
        PTL(2) = -.880239153726985902;
        PTL(3) = -.78151400389680141;
        PTL(4) = -.657671159216690766;
        PTL(5) = -.512690537086476968;
        PTL(6) = -.351231763453876315;
        PTL(7) = -.178484181495847856;
        PTL(8) = .0;

        PTW(0) = .024148302868547932;
        PTW(1) = .055459529373987201;
        PTW(2) = .085036148317179181;
        PTW(3) = .11188384719340397;
        PTW(4) = .135136368468525473;
        PTW(5) = .15404576107681029;
        PTW(6) = .168004102156450045;
        PTW(7) = .176562705366992646;
        PTW(8) = .179446470356206526;

        break;
    }
    case 18: {
        PTL(0) = -.991565168420930947;
        PTL(1) = -.955823949571397755;
        PTL(2) = -.892602466497555739;
        PTL(3) = -.803704958972523116;
        PTL(4) = -.69168704306035321;
        PTL(5) = -.559770831073947535;
        PTL(6) = -.411751161462842646;
        PTL(7) = -.25188622569150551;
        PTL(8) = -.084775013041735301;

        PTW(0) = .02161601352648331;
        PTW(1) = .0497145488949698;
        PTW(2) = .076425730254889057;
        PTW(3) = .100942044106287166;
        PTW(4) = .12255520671147846;
        PTW(5) = .140642914670650651;
        PTW(6) = .154684675126265245;
        PTW(7) = .164276483745832723;
        PTW(8) = .169142382963143592;

        break;
    }
    case 19: {
        PTL(0) = -.992406843843584403;
        PTL(1) = -.960208152134830031;
        PTL(2) = -.903155903614817902;
        PTL(3) = -.82271465653714283;
        PTL(4) = -.720966177335229379;
        PTL(5) = -.600545304661681024;
        PTL(6) = -.464570741375960946;
        PTL(7) = -.316564099963629832;
        PTL(8) = -.160358645640225376;
        PTL(9) = .0;

        PTW(0) = .019461788229726477;
        PTW(1) = .0448142267656996;
        PTW(2) = .069044542737641227;
        PTW(3) = .091490021622449999;
        PTW(4) = .111566645547333995;
        PTW(5) = .128753962539336228;
        PTW(6) = .142606702173606612;
        PTW(7) = .152766042065859667;
        PTW(8) = .15896884339395435;
        PTW(9) = .161054449848783696;

        break;
    }
    case 20: {
        PTL(0) = -.993128599185094925;
        PTL(1) = -.963971927277913791;
        PTL(2) = -.912234428251325906;
        PTL(3) = -.839116971822218823;
        PTL(4) = -.746331906460150793;
        PTL(5) = -.636053680726515026;
        PTL(6) = -.510867001950827098;
        PTL(7) = -.373706088715419561;
        PTL(8) = -.227785851141645078;
        PTL(9) = -.0765265211334973338;

        PTW(0) = .017614007139152118;
        PTW(1) = .040601429800386941;
        PTW(2) = .062672048334109064;
        PTW(3) = .083276741576704749;
        PTW(4) = .10193011981724044;
        PTW(5) = .118194531961518417;
        PTW(6) = .131688638449176627;
        PTW(7) = .142096109318382051;
        PTW(8) = .149172986472603747;
        PTW(9) = .152753387130725851;

        break;
    }
    default:
        throw std::invalid_argument("not supported");
    }

    for(auto I = 0u, J = order - 1; I < order / 2; ++I, --J) {
        PTL(J) = -PTL(I);
        PTW(J) = PTW(I);
    }

    return {std::move(PTL), std::move(PTW)};
}

template<> std::pair<arma::vec, arma::vec> generate_seeds<IntegrationType::HERMITE>(const unsigned order) {
    arma::vec PTL(order, arma::fill::none), PTW(order, arma::fill::none);

    switch(order) {
    case 2: {
        PTL(0) = -0.7071067811;

        PTW(0) = 0.8862269254;

        break;
    }
    case 3: {
        PTL(0) = -1.2247448713;
        PTL(1) = 0.;

        PTW(0) = 0.2954089751;
        PTW(1) = 1.1816359006;

        break;
    }
    case 4: {
        PTL(0) = -1.6506801238;
        PTL(1) = -0.5246476232;

        PTW(0) = 0.0813128354;
        PTW(1) = 0.8049140900;

        break;
    }
    case 5: {
        PTL(0) = -2.0201828704;
        PTL(1) = -0.9585724646;
        PTL(2) = 0.;

        PTW(0) = 0.0199532420;
        PTW(1) = 0.3936193231;
        PTW(2) = 0.9453087204;

        break;
    }
    case 6: {
        PTL(0) = -2.3506049736;
        PTL(1) = -1.3358490740;
        PTL(2) = -0.4360774119;

        PTW(0) = 0.0045300099;
        PTW(1) = 0.1570673203;
        PTW(2) = 0.7246295952;

        break;
    }
    default:
        throw std::invalid_argument("not supported");
    }

    for(auto I = 0u, J = order - 1; I < order / 2; ++I, --J) {
        PTL(J) = -PTL(I);
        PTW(J) = PTW(I);
    }

    return {std::move(PTL), std::move(PTW)};
}

template<> std::pair<arma::vec, arma::vec> generate_seeds<IntegrationType::CHEBYSHEV>(const unsigned order) {
    arma::vec PTL(order, arma::fill::none), PTW(order, arma::fill::none);

    switch(order) {
    case 2: {
        PTL(0) = -0.5773502691;

        PTW(0) = 1.;

        break;
    }
    case 3: {
        PTL(0) = -0.7071067811;
        PTL(1) = 0.;

        PTW(0) = 2. / 3.;
        PTW(1) = 2. / 3.;

        break;
    }
    case 4: {
        PTL(0) = -0.7946544722;
        PTL(1) = -0.1875924740;

        PTW(0) = .5;
        PTW(1) = .5;

        break;
    }
    case 5: {
        PTL(0) = -0.8324974870;
        PTL(1) = -0.3745414095;
        PTL(2) = 0.;

        PTW(0) = .4;
        PTW(1) = .4;
        PTW(2) = .4;

        break;
    }
    case 6: {
        PTL(0) = -0.8662468181;
        PTL(1) = -0.4225186537;
        PTL(2) = -0.2666354015;

        PTW(0) = 1. / 3.;
        PTW(1) = 1. / 3.;
        PTW(2) = 1. / 3.;

        break;
    }
    case 7: {
        PTL(0) = -0.8838617007;
        PTL(1) = -0.5296567752;
        PTL(2) = -0.3239118105;
        PTL(3) = 0.;

        PTW(0) = 2. / 7.;
        PTW(1) = 2. / 7.;
        PTW(2) = 2. / 7.;
        PTW(3) = 2. / 7.;

        break;
    }
    default:
        throw std::invalid_argument("not supported");
    }

    for(auto I = 0u, J = order - 1; I < order / 2; ++I, --J) {
        PTL(J) = -PTL(I);
        PTW(J) = PTW(I);
    }

    return {std::move(PTL), std::move(PTW)};
}

template<> std::pair<arma::vec, arma::vec> generate_seeds<IntegrationType::LOBATTO>(const unsigned order) {
    arma::vec PTL(order, arma::fill::none), PTW(order, arma::fill::none);

    switch(order) {
    case 3: {
        PTL(0) = -1.;
        PTL(1) = .0;

        PTW(0) = 1. / 3.;
        PTW(1) = 4. / 3.;

        break;
    }
    case 4: {
        PTL(0) = -1.;
        PTL(1) = -sqrt(.2);

        PTW(0) = 1. / 6.;
        PTW(1) = 5. / 6.;

        break;
    }
    case 5: {
        PTL(0) = -1.;
        PTL(1) = -sqrt(3. / 7.);
        PTL(2) = .0;

        PTW(0) = .1;
        PTW(1) = 49. / 90.;
        PTW(2) = 64. / 90.;

        break;
    }
    case 6: {
        const auto TMPA = 2. * sqrt(7.);
        constexpr auto TMPB = 14. / 30.;
        const auto TMPC = sqrt(7.) / 30.;

        PTL(0) = -1.;
        PTL(1) = -sqrt((7. + TMPA) / 21.);
        PTL(2) = -sqrt((7. - TMPA) / 21.);

        PTW(0) = 1. / 15.;
        PTW(1) = TMPB - TMPC;
        PTW(2) = TMPB + TMPC;

        break;
    }
    case 7: {
        PTL(0) = -1.;
        PTL(1) = -.83022389627856693;
        PTL(2) = -.468848793470714214;
        PTL(3) = .0;

        PTW(0) = .0476190476190476191;
        PTW(1) = .276826047361565948;
        PTW(2) = .431745381209862623;
        PTW(3) = .48761904761904762;

        break;
    }
    case 8: {
        PTL(0) = -1.;
        PTL(1) = -.871740148509606615;
        PTL(2) = -.591700181433142302;
        PTL(3) = -.209299217902478869;

        PTW(0) = .0357142857142857143;
        PTW(1) = .210704227143506039;
        PTW(2) = .34112269248350436;
        PTW(3) = .412458794658703882;

        break;
    }
    case 9: {
        PTL(0) = -1.;
        PTL(1) = -.899757995411460157;
        PTL(2) = -.677186279510737753;
        PTL(3) = -.363117463826178159;
        PTL(4) = .0;

        PTW(0) = .0277777777777777778;
        PTW(1) = .165495361560805525;
        PTW(2) = .274538712500161735;
        PTW(3) = .346428510973046345;
        PTW(4) = .371519274376417234;

        break;
    }
    case 10: {
        PTL(0) = -1.;
        PTL(1) = -.919533908166458814;
        PTL(2) = -.738773865105505075;
        PTL(3) = -.477924949810444496;
        PTL(4) = -.165278957666387025;

        PTW(0) = .0222222222222222222;
        PTW(1) = .13330599085107011;
        PTW(2) = .22488934206312645;
        PTW(3) = .292042683679683758;
        PTW(4) = .327539761183897457;

        break;
    }
    case 11: {
        PTL(0) = -1.;
        PTL(1) = -.934001430408059134;
        PTL(2) = -.78448347366314442;
        PTL(3) = -.565235326996205007;
        PTL(4) = -.295758135586939391;
        PTL(5) = .0;

        PTW(0) = .0181818181818181818;
        PTW(1) = .109612273266994865;
        PTW(2) = .187169881780305204;
        PTW(3) = .248048104264028314;
        PTW(4) = .286879124779008089;
        PTW(5) = .30021759545569069;

        break;
    }
    case 12: {
        PTL(0) = -1.;
        PTL(1) = -.944899272222882223;
        PTL(2) = -.819279321644006678;
        PTL(3) = -.632876153031860678;
        PTL(4) = -.399530940965348932;
        PTL(5) = -.136552932854927555;

        PTW(0) = .0151515151515151515;
        PTW(1) = .09168451741319613;
        PTW(2) = .15797470556437012;
        PTW(3) = .212508417761021145;
        PTW(4) = .25127560319920128;
        PTW(5) = .271405240910696177;

        break;
    }
    case 13: {
        PTL(0) = -1.;
        PTL(1) = -.953309846642163912;
        PTL(2) = -.84634756465187232;
        PTL(3) = -.686188469081757426;
        PTL(4) = -.482909821091336202;
        PTL(5) = -.249286930106239993;
        PTL(6) = .0;

        PTW(0) = .0128205128205128205;
        PTW(1) = .07780168674681893;
        PTW(2) = .134981926689608349;
        PTW(3) = .18364686520355009;
        PTW(4) = .220767793566110086;
        PTW(5) = .24401579030667636;
        PTW(6) = .251930849333446736;

        break;
    }
    case 14: {
        PTL(0) = -1.;
        PTL(1) = -.959935045267260901;
        PTL(2) = -.867801053830347251;
        PTL(3) = -.728868599091326141;
        PTL(4) = -.550639402928647055;
        PTL(5) = -.342724013342712845;
        PTL(6) = -.116331868883703868;

        PTW(0) = .010989010989010989;
        PTW(1) = .06683728449768128;
        PTW(2) = .116586655898711652;
        PTW(3) = .160021851762952142;
        PTW(4) = .19482614937341612;
        PTW(5) = .219126253009770755;
        PTW(6) = .231612794468457059;

        break;
    }
    case 15: {
        PTL(0) = -1.;
        PTL(1) = -.965245926503838573;
        PTL(2) = -.885082044222976299;
        PTL(3) = -.763519689951815201;
        PTL(4) = -.606253205469845711;
        PTL(5) = -.420638054713672481;
        PTL(6) = -.215353955363794238;
        PTL(7) = .0;

        PTW(0) = .00952380952380952381;
        PTW(1) = .05802989302860125;
        PTW(2) = .101660070325718068;
        PTW(3) = .14051169980242811;
        PTW(4) = .172789647253600949;
        PTW(5) = .19698723596461336;
        PTW(6) = .21197358592682092;
        PTW(7) = .21704811634881565;

        break;
    }
    case 16: {
        PTL(0) = -1.;
        PTL(1) = -.969568046270217933;
        PTL(2) = -.899200533093472093;
        PTL(3) = -.792008291861815064;
        PTL(4) = -.65238870288249309;
        PTL(5) = -.486059421887137612;
        PTL(6) = -.299830468900763208;
        PTL(7) = -.101326273521949448;

        PTW(0) = .00833333333333333333;
        PTW(1) = .05085036100591991;
        PTW(2) = .089393697325930801;
        PTW(3) = .124255382132514098;
        PTW(4) = .154026980807164281;
        PTW(5) = .177491913391704125;
        PTW(6) = .19369002382520358;
        PTW(7) = .201958308178229872;

        break;
    }
    case 17: {
        PTL(0) = -1.;
        PTL(1) = -.973132176631418314;
        PTL(2) = -.910879995915573596;
        PTL(3) = -.815696251221770307;
        PTL(4) = -.691028980627684705;
        PTL(5) = -.541385399330101539;
        PTL(6) = -.372174433565477042;
        PTL(7) = -.189511973518317388;
        PTL(8) = .0;

        PTW(0) = .00735294117647058824;
        PTW(1) = .04492194054325421;
        PTW(2) = .079198270503687119;
        PTW(3) = .110592909007028161;
        PTW(4) = .137987746201926559;
        PTW(5) = .16039466199762154;
        PTW(6) = .17700425351565787;
        PTW(7) = .18721633967761924;
        PTW(8) = .190661874753469433;

        break;
    }
    case 18: {
        PTL(0) = -1.;
        PTL(1) = -.976105557412198543;
        PTL(2) = -.920649185347533874;
        PTL(3) = -.835593535218090214;
        PTL(4) = -.723679329283242681;
        PTL(5) = -.588504834318661761;
        PTL(6) = -.434415036912123975;
        PTL(7) = -.26636265287828098;
        PTL(8) = -.089749093484652111;

        PTW(0) = .00653594771241830065;
        PTW(1) = .039970628810914066;
        PTW(2) = .070637166885633665;
        PTW(3) = .0990162717175028;
        PTW(4) = .1242105331329671;
        PTW(5) = .145411961573802268;
        PTW(6) = .161939517237602489;
        PTW(7) = .173262109489456226;
        PTW(8) = .17901586343970308;

        break;
    }
    case 19: {
        PTL(0) = -1.;
        PTL(1) = -.978611766222080095;
        PTL(2) = -.928901528152586244;
        PTL(3) = -.852460577796646093;
        PTL(4) = -.751494202552613014;
        PTL(5) = -.628908137265220498;
        PTL(6) = -.488229285680713503;
        PTL(7) = -.33350484782449861;
        PTL(8) = -.169186023409281571;
        PTL(9) = .0;

        PTW(0) = .00584795321637426901;
        PTW(1) = .035793365186176477;
        PTW(2) = .063381891762629737;
        PTW(3) = .0891317570992070845;
        PTW(4) = .11231534147730504;
        PTW(5) = .132267280448750777;
        PTW(6) = .14841394259593889;
        PTW(7) = .160290924044061242;
        PTW(8) = .167556584527142867;
        PTW(9) = .170001919284827235;

        break;
    }
    case 20: {
        PTL(0) = -1.;
        PTL(1) = -.980743704893914172;
        PTL(2) = -.935934498812665436;
        PTL(3) = -.866877978089950141;
        PTL(4) = -.77536826095205587;
        PTL(5) = -.66377640229031129;
        PTL(6) = -.53499286403188626;
        PTL(7) = -.392353183713909299;
        PTL(8) = -.239551705922986495;
        PTL(9) = -.080545937238821838;

        PTW(0) = .00526315789473684211;
        PTW(1) = .032237123188488941;
        PTW(2) = .057181802127566826;
        PTW(3) = .080631763996119603;
        PTW(4) = .10199149969945082;
        PTW(5) = .12070922762867473;
        PTW(6) = .136300482358724185;
        PTW(7) = .148361554070916826;
        PTW(8) = .156580102647475487;
        PTW(9) = .16074328638784575;

        break;
    }
    default:
        throw std::invalid_argument("not supported");
    }

    for(auto I = 0u, J = order - 1; I < order / 2; ++I, --J) {
        PTL(J) = -PTL(I);
        PTW(J) = PTW(I);
    }

    return {std::move(PTL), std::move(PTW)};
}

template<> std::pair<arma::vec, arma::vec> generate_seeds<IntegrationType::RADAU>(const unsigned order) {
    arma::vec PTL(order, arma::fill::none), PTW(order, arma::fill::none);

    switch(order) {
    case 2: {
        PTL(0) = -1.;
        PTL(1) = 1. / 3.;

        PTW(0) = .5;
        PTW(1) = 1.5;

        break;
    }
    case 3: {
        PTL(0) = -1.;
        PTL(1) = -.2898979485;
        PTL(2) = .6898979485;

        PTW(0) = 2. / 9.;
        PTW(1) = 1.0249716523;
        PTW(2) = .7528061254;

        break;
    }
    case 4: {
        PTL(0) = -1.;
        PTL(1) = -.5753189235;
        PTL(2) = .1810662711;
        PTL(3) = .8228240809;

        PTW(0) = .125;
        PTW(1) = .6576886399;
        PTW(2) = .7763869376;
        PTW(3) = .4409244223;

        break;
    }
    case 5: {
        PTL(0) = -1.;
        PTL(1) = -.7204802713;
        PTL(2) = -.1671808647;
        PTL(3) = .44631398727;
        PTL(4) = .8857916077;

        PTW(0) = .08;
        PTW(1) = .4462078021;
        PTW(2) = .6236530459;
        PTW(3) = .5627120302;
        PTW(4) = .2874271215;

        break;
    }
    case 6: {
        PTL(0) = -1.;
        PTL(1) = -.8029298284;
        PTL(2) = -.3909285467;
        PTL(3) = .1240503795;
        PTL(4) = .6039731642;
        PTL(5) = .9203802858;

        PTW(0) = 1. / 18.;
        PTW(1) = .3196407532;
        PTW(2) = .4853871884;
        PTW(3) = .5209267831;
        PTW(4) = .4169013343;
        PTW(5) = .2015883852;

        break;
    }
    case 7: {
        PTL(0) = -1.;
        PTL(1) = -.8538913426;
        PTL(2) = -.5384677240;
        PTL(3) = -.1173430375;
        PTL(4) = .3260306194;
        PTL(5) = .7038428006;
        PTL(6) = .9413671456;

        PTW(0) = 2. / 49.;
        PTW(1) = .2392274892;
        PTW(2) = .3809498736;
        PTW(3) = .4471098290;
        PTW(4) = .4247037790;
        PTW(5) = .3182042314;
        PTW(6) = .1489884711;

        break;
    }
    case 8: {
        PTL(0) = -1.;
        PTL(1) = -.8874748789;
        PTL(2) = -.6395186165;
        PTL(3) = -.2947505657;
        PTL(4) = .0943072526;
        PTL(5) = .4684203544;
        PTL(6) = .7706418936;
        PTL(7) = .9550412271;

        PTW(0) = 1. / 32.;
        PTW(1) = .1853581548;
        PTW(2) = .3041306206;
        PTW(3) = .3765175453;
        PTW(4) = .3915721674;
        PTW(5) = .3470147956;
        PTW(6) = .2496479013;
        PTW(7) = .1145088147;

        break;
    }
    case 9: {
        PTL(0) = -1.;
        PTL(1) = -.9107320894;
        PTL(2) = -.7112674859;
        PTL(3) = -.4263504857;
        PTL(4) = -.0903733696;
        PTL(5) = .2561356708;
        PTL(6) = .5713830412;
        PTL(7) = .8173527842;
        PTL(8) = .9644401697;

        PTW(0) = 2. / 81.;
        PTW(1) = .1476540190;
        PTW(2) = .2471893782;
        PTW(3) = .3168437756;
        PTW(4) = .3482730027;
        PTW(5) = .3376939669;
        PTW(6) = .2863866963;
        PTW(7) = .2005532980;
        PTW(8) = .0907145049;

        break;
    }
    case 10: {
        PTL(0) = -1.;
        PTL(1) = -.9274843742;
        PTL(2) = -.7638420424;
        PTL(3) = -.5256460303;
        PTL(4) = -.2362344693;
        PTL(5) = .0760591978;
        PTL(6) = .3806648401;
        PTL(7) = .6477666876;
        PTL(8) = .8512252205;
        PTL(9) = .9711751807;

        PTW(0) = .02;
        PTW(1) = .1202966705;
        PTW(2) = .2042701318;
        PTW(3) = .2681948378;
        PTW(4) = .3058592877;
        PTW(5) = .3135824572;
        PTW(6) = .2906101648;
        PTW(7) = .2391934317;
        PTW(8) = .1643760127;
        PTW(9) = .0736170054;

        break;
    }
    default:
        throw std::invalid_argument("not supported");
    }

    return {std::move(PTL), std::move(PTW)};
}

template<> std::pair<arma::vec, arma::vec> generate_seeds<IntegrationType::LAGUERRE>(const unsigned order) {
    arma::vec PTL(order, arma::fill::none), PTW(order, arma::fill::none);

    switch(order) {
    case 2: {
        PTL(0) = 0.5857864376;
        PTL(1) = 3.4142135623;

        PTW(0) = 0.8535533905;
        PTW(1) = 0.1464466094;

        break;
    }
    case 3: {
        PTL(0) = 0.4157745567;
        PTL(1) = 2.2942803602;
        PTL(2) = 6.2899450829;

        PTW(0) = 0.7110930099;
        PTW(1) = 0.2785177335;
        PTW(2) = 0.0103892565;

        break;
    }
    case 4: {
        PTL(0) = 0.3225476896;
        PTL(1) = 1.7457611011;
        PTL(2) = 4.5366202969;
        PTL(3) = 9.3950709123;

        PTW(0) = 0.6031541043;
        PTW(1) = 0.3574186924;
        PTW(2) = 0.0388879085;
        PTW(3) = 0.0005392947;

        break;
    }
    case 5: {
        PTL(0) = 0.2635603197;
        PTL(1) = 1.4134030591;
        PTL(2) = 3.5964257710;
        PTL(3) = 7.0858100058;
        PTL(4) = 12.6408008442;

        PTW(0) = 0.5217556105;
        PTW(1) = 0.3986668110;
        PTW(2) = 0.0759424496;
        PTW(3) = 0.0036117586;
        PTW(4) = 0.0000233699;

        break;
    }
    default:
        throw std::invalid_argument("not supported");
    }

    return {std::move(PTL), std::move(PTW)};
}

IntegrationPlan::IntegrationPlan(const unsigned D, const unsigned O, const IntegrationType T) {
    switch(T) {
    case IntegrationType::GAUSS:
        generate<IntegrationType::GAUSS>(D, O);
        break;
    case IntegrationType::HERMITE:
        generate<IntegrationType::HERMITE>(D, O);
        break;
    case IntegrationType::CHEBYSHEV:
        generate<IntegrationType::CHEBYSHEV>(D, O);
        break;
    case IntegrationType::LOBATTO:
        generate<IntegrationType::LOBATTO>(D, O);
        break;
    case IntegrationType::RADAU:
        generate<IntegrationType::RADAU>(D, O);
        break;
    case IntegrationType::LAGUERRE:
        generate<IntegrationType::LAGUERRE>(D, O);
        break;
    case IntegrationType::IRONS:
        generate<IntegrationType::IRONS>(D, O);
        break;
    case IntegrationType::TRIANGLE:
        generate<IntegrationType::TRIANGLE>(D, O);
        break;
    }
}

const arma::mat& IntegrationPlan::get_data() const { return int_pts; }

double IntegrationPlan::operator()(const unsigned i, const unsigned j) const { return int_pts(i, j); }

void IntegrationPlan::print() const {
    for(unsigned i = 0; i < n_rows; ++i) {
        printf("Node %u\t", i + 1);
        for(unsigned j = 0; j < n_cols - 1; ++j) printf("%+.6E\t", int_pts(i, j));
        printf("Weight\t%+.6E\n", int_pts(i, n_cols - 1));
    }
    printf("\n");
}

template<> void IntegrationPlan::generate<IntegrationType::GAUSS>(const unsigned D, const unsigned O) {
    int_pts = generate_points(D, generate_seeds<IntegrationType::GAUSS>(std::min(std::max(O, 1u), 20u)));
    arma::access::rw(n_rows) = static_cast<unsigned>(int_pts.n_rows);
    arma::access::rw(n_cols) = static_cast<unsigned>(int_pts.n_cols);
}

template<> void IntegrationPlan::generate<IntegrationType::HERMITE>(const unsigned D, const unsigned O) {
    int_pts = generate_points(D, generate_seeds<IntegrationType::HERMITE>(std::min(std::max(O, 2u), 6u)));
    arma::access::rw(n_rows) = static_cast<unsigned>(int_pts.n_rows);
    arma::access::rw(n_cols) = static_cast<unsigned>(int_pts.n_cols);
}

template<> void IntegrationPlan::generate<IntegrationType::CHEBYSHEV>(const unsigned D, const unsigned O) {
    int_pts = generate_points(D, generate_seeds<IntegrationType::CHEBYSHEV>(std::min(std::max(O, 2u), 7u)));
    arma::access::rw(n_rows) = static_cast<unsigned>(int_pts.n_rows);
    arma::access::rw(n_cols) = static_cast<unsigned>(int_pts.n_cols);
}

template<> void IntegrationPlan::generate<IntegrationType::LOBATTO>(const unsigned D, const unsigned O) {
    int_pts = generate_points(D, generate_seeds<IntegrationType::LOBATTO>(std::min(std::max(O, 3u), 20u)));
    arma::access::rw(n_rows) = static_cast<unsigned>(int_pts.n_rows);
    arma::access::rw(n_cols) = static_cast<unsigned>(int_pts.n_cols);
}

template<> void IntegrationPlan::generate<IntegrationType::RADAU>(const unsigned D, const unsigned O) {
    int_pts = generate_points(D, generate_seeds<IntegrationType::RADAU>(std::min(std::max(O, 2u), 10u)));
    arma::access::rw(n_rows) = static_cast<unsigned>(int_pts.n_rows);
    arma::access::rw(n_cols) = static_cast<unsigned>(int_pts.n_cols);
}

template<> void IntegrationPlan::generate<IntegrationType::LAGUERRE>(const unsigned D, const unsigned O) {
    int_pts = generate_points(D, generate_seeds<IntegrationType::LAGUERRE>(std::min(std::max(O, 2u), 5u)));
    arma::access::rw(n_rows) = static_cast<unsigned>(int_pts.n_rows);
    arma::access::rw(n_cols) = static_cast<unsigned>(int_pts.n_cols);
}

template<> void IntegrationPlan::generate<IntegrationType::IRONS>(const unsigned D, const unsigned O) {
    arma::access::rw(n_cols) = D + 1;

    if(3 == D && 2 == O) {
        arma::access::rw(n_rows) = 6;

        int_pts.set_size(n_rows, n_cols);

        constexpr auto WB = 4. / 3.;
        for(auto I = 0; I < 6; ++I) int_pts(I, 3) = WB;
        int_pts(0, 0) = int_pts(2, 1) = int_pts(4, 2) = -(int_pts(1, 0) = int_pts(3, 1) = int_pts(5, 2) = 1.);
        int_pts(0, 1) = int_pts(0, 2) = int_pts(1, 1) = int_pts(1, 2) = int_pts(2, 0) = int_pts(2, 2) = int_pts(3, 0) = int_pts(3, 2) = int_pts(4, 0) = int_pts(4, 1) = int_pts(5, 0) = int_pts(5, 1) = 0.;

        return;
    }

    if(3 == D && 3 == O) {
        arma::access::rw(n_rows) = 14;

        int_pts.set_size(n_rows, n_cols);

        constexpr auto WB = .886426593;
        constexpr auto WC = .335180055;
        constexpr auto LB = .795822426;
        constexpr auto LC = .758786911;
        for(auto I = 0; I < 6; ++I) int_pts(I, 3) = WB;
        for(auto I = 6; I < 14; ++I) int_pts(I, 3) = WC;
        int_pts(0, 0) = int_pts(2, 1) = int_pts(4, 2) = -(int_pts(1, 0) = int_pts(3, 1) = int_pts(5, 2) = LB);
        int_pts(0, 1) = int_pts(0, 2) = int_pts(1, 1) = int_pts(1, 2) = int_pts(2, 0) = int_pts(2, 2) = int_pts(3, 0) = int_pts(3, 2) = int_pts(4, 0) = int_pts(4, 1) = int_pts(5, 0) = int_pts(5, 1) = 0.;
        int_pts(6, 0) = int_pts(6, 1) = int_pts(6, 2) = int_pts(7, 1) = int_pts(7, 2) = int_pts(8, 2) = int_pts(9, 0) = int_pts(9, 2) = int_pts(10, 0) = int_pts(10, 1) = int_pts(11, 1) = int_pts(13, 0) = -LC;
        int_pts(7, 0) = int_pts(8, 0) = int_pts(8, 1) = int_pts(9, 1) = int_pts(10, 2) = int_pts(11, 0) = int_pts(11, 2) = int_pts(12, 0) = int_pts(12, 1) = int_pts(12, 2) = int_pts(13, 1) = int_pts(13, 2) = LC;

        return;
    }

    if(2 == D && 2 == O) {
        arma::access::rw(n_rows) = 5;

        int_pts.set_size(n_rows, n_cols);

        constexpr auto WB = 2. / 3.;
        for(auto I = 0; I < 5; ++I) {
            for(auto J = 0; J < 2; ++J) int_pts(I, J) = 0.;
            int_pts(I, 2) = WB;
        }
        int_pts(0, 2) *= 2.;
        int_pts(1, 0) = int_pts(2, 1) = -1.;
        int_pts(3, 0) = int_pts(4, 1) = 1.;

        return;
    }

    throw std::invalid_argument("not supported");
}

template<> void IntegrationPlan::generate<IntegrationType::TRIANGLE>(const unsigned D, const unsigned O) {
    if(3 == D) throw std::invalid_argument("not supported");

    arma::access::rw(n_cols) = 4;

    if(1 == O) {
        arma::access::rw(n_rows) = 1;

        int_pts.set_size(n_rows, n_cols);

        int_pts(0, 0) = int_pts(0, 1) = int_pts(0, 2) = 1. / 3.;
        int_pts(0, 3) = 1.;

        return;
    }

    if(2 == O) {
        arma::access::rw(n_rows) = 3;

        int_pts.set_size(n_rows, n_cols);

        int_pts(0, 0) = int_pts(0, 1) = int_pts(1, 1) = int_pts(1, 2) = int_pts(2, 0) = int_pts(2, 2) = .5;
        int_pts(0, 2) = int_pts(1, 0) = int_pts(2, 1) = 0.;
        int_pts(0, 3) = int_pts(1, 3) = int_pts(2, 3) = 1. / 3.;

        return;
    }

    if(3 == O) {
        arma::access::rw(n_rows) = 4;

        int_pts.set_size(n_rows, n_cols);

        int_pts(0, 0) = int_pts(0, 1) = int_pts(0, 2) = 1. / 3.;
        int_pts(0, 3) = -27. / 48.;

        int_pts(1, 1) = int_pts(1, 2) = int_pts(2, 0) = int_pts(2, 2) = int_pts(3, 0) = int_pts(3, 1) = .2;
        int_pts(1, 0) = int_pts(2, 1) = int_pts(3, 2) = .6;
        int_pts(1, 3) = int_pts(2, 3) = int_pts(3, 3) = 25. / 48.;

        return;
    }

    if(5 == O) {
        arma::access::rw(n_rows) = 7;

        int_pts.set_size(n_rows, n_cols);

        int_pts(0, 0) = int_pts(0, 1) = int_pts(0, 2) = 1. / 3.;
        int_pts(0, 3) = .225;

        int_pts(1, 1) = int_pts(1, 2) = int_pts(2, 0) = int_pts(2, 2) = int_pts(3, 0) = int_pts(3, 1) = .4701420641;
        int_pts(1, 0) = int_pts(2, 1) = int_pts(3, 2) = .0597158717;
        int_pts(1, 3) = int_pts(2, 3) = int_pts(3, 3) = .1323941527;

        int_pts(4, 1) = int_pts(4, 2) = int_pts(5, 0) = int_pts(5, 2) = int_pts(6, 0) = int_pts(6, 1) = .1012865073;
        int_pts(4, 0) = int_pts(5, 1) = int_pts(6, 2) = .7974269853;
        int_pts(4, 3) = int_pts(5, 3) = int_pts(6, 3) = .1259391805;

        return;
    }

    throw std::invalid_argument("not supported");
}
