/*
 * Test the AABBTree functionality.
 */

#include <random>
#include <utility>
#include <iostream>
#include <aabbtree.hpp>

namespace flottekarte {

static std::vector<segment_t> test_data()
{
    std::vector<segment_t> data;
    data.emplace_back(xy_t(0.100734613189731703, 0.46137984678539673),
                      xy_t(0.0417424656059637578, 0.543239871991812229));
    data.emplace_back(xy_t(0.61504936908344443, 0.709263361612838428),
                      xy_t(0.626213515674343135, 0.00790235077145155784));
    data.emplace_back(xy_t(0.592909810307182039, 0.134776440318631452),
                      xy_t(0.341908154623234928, 0.652257240116251302));
    data.emplace_back(xy_t(0.680118719264586558, 0.133870361392252113),
                      xy_t(0.528694827668636846, 0.329930575422497763));
    data.emplace_back(xy_t(0.590126016170673062, 0.956435155469485321),
                      xy_t(0.933529493327886262, 0.489735247443793431));
    data.emplace_back(xy_t(0.536378227229399762, 0.617245529289275896),
                      xy_t(0.866586264061118716, 0.307129141810447426));
    data.emplace_back(xy_t(0.790238834243362387, 0.590169080083150277),
                      xy_t(0.390602334247182215, 0.349787900973318189));
    data.emplace_back(xy_t(0.81411161282537825, 0.266529492669613122),
                      xy_t(0.743837964713740418, 0.372136412322686072));
    data.emplace_back(xy_t(0.292699033093345051, 0.934418276380332702),
                      xy_t(0.784378549208133324, 0.497274806615162379));
    data.emplace_back(xy_t(0.343697197795442999, 0.65070754700971456),
                      xy_t(0.32019587886077111, 0.247121443453660022));
    data.emplace_back(xy_t(0.717884614234240215, 0.605211105448415609),
                      xy_t(0.339401918898149801, 0.556092984802551094));
    data.emplace_back(xy_t(0.824011780822617057, 0.498627233125908509),
                      xy_t(0.73603784841780262, 0.093194160458210773));
    data.emplace_back(xy_t(0.983577368035390309, 0.548428138637450102),
                      xy_t(0.522980314288423065, 0.933665397126846486));
    data.emplace_back(xy_t(0.289389561575763399, 0.243399004163558502),
                      xy_t(0.608311651715353841, 0.891933716725878845));
    data.emplace_back(xy_t(0.361784247860436892, 0.830943101205916967),
                      xy_t(0.321211175427722595, 0.451427588087141396));
    data.emplace_back(xy_t(0.374644653297708796, 0.371096458025590548),
                      xy_t(0.484125571396318732, 0.495697268982999095));
    data.emplace_back(xy_t(0.160267337478929461, 0.108063542728441916),
                      xy_t(0.547672883410344591, 0.806656832464973661));
    data.emplace_back(xy_t(0.0764033164226872036, 0.750848925126123445),
                      xy_t(0.916839711441964855, 0.141924494363830811));
    data.emplace_back(xy_t(0.204350389465422527, 0.68896476048738331),
                      xy_t(0.146773187243761732, 0.884290907010028215));
    data.emplace_back(xy_t(0.949542244840127947, 0.644135199915178869),
                      xy_t(0.281170392512216094, 0.19549350762180312));
    data.emplace_back(xy_t(0.707032605675405823, 0.560413916762580011),
                      xy_t(0.631893761208098925, 0.924924366551614585));
    data.emplace_back(xy_t(0.56224215625614371, 0.73229320212639526),
                      xy_t(0.691476373539519518, 0.26249607308860462));
    data.emplace_back(xy_t(0.00129594475605649198, 0.41834769748046341),
                      xy_t(0.966542574886245953, 0.0585905242343896862));
    data.emplace_back(xy_t(0.763203585837627263, 0.0136939880611458929),
                      xy_t(0.948795729957973677, 0.426695982755321002));
    data.emplace_back(xy_t(0.761755817745130015, 0.909612246243156175),
                      xy_t(0.671208930983029139, 0.957954138088308205));
    data.emplace_back(xy_t(0.303993918455004619, 0.880726220084155909),
                      xy_t(0.350701143736570076, 0.213887533758673365));
    data.emplace_back(xy_t(0.150603491481951013, 0.149341821318858509),
                      xy_t(0.20437592064542362, 0.853526640186479058));
    data.emplace_back(xy_t(0.192741503120031826, 0.267023009692266056),
                      xy_t(0.744238700638989892, 0.738403196292783393));
    data.emplace_back(xy_t(0.517675041090769228, 0.338649658277970056),
                      xy_t(0.109232785315109604, 0.0800899162733100384));
    data.emplace_back(xy_t(0.0354074961751647682, 0.509608018799408979),
                      xy_t(0.534375093007446234, 0.10656650307474054));

    return data;
}

/*
 * The function that generated that test data, not used here:
 */
static void generate_test_data()
{

    std::default_random_engine rng(8299832);
    std::uniform_real_distribution dist(0.0, 1.0);
    std::vector<segment_t> segments;
    std::cout << std::setprecision(18);
    for (size_t i=0; i<30; ++i)
    {
        segments.emplace_back(xy_t(dist(rng), dist(rng)),
                              xy_t(dist(rng), dist(rng)));
        std::cout << "    data.emplace_back(xy_t(" << segments.back().p0.x << ", "
            << segments.back().p0.y << "),\n"
               "                      xy_t("
            << segments.back().p1.x << ", "
            << segments.back().p1.y << "));\n";
    }
    std::cout << std::flush;
}



/*
 * A brute force intersection algorithm.
 */
template<typename geometry_t>
bool brute_force_intersection(
    const geometry_t& geom,
    const std::vector<segment_t>& segments
)
{
    for (size_t i=0; i<segments.size(); ++i){
        if (intersects(geom, segments[i])){
            std::cout << "intersection at segment " << i << "\n";
            return true;
        }
    }
    return false;
}

#pragma GCC push_options
#pragma GCC optimize ("O0")
void test_aabbtree()
{
//    generate_test_data();

    AABBTree tree;

    /* The intersection test geometries: */
    const segment_t the_segment(xy_t(0.07, 0.07),xy_t(0.2, 0.27));
    const circle_t the_circle(xy_t(0.4, 0.4), 0.02);
    const circle_t circle2(xy_t(0.4, 0.4), 0.05);

    /* Load the randomly generated segment cloud: */
    std::vector<segment_t> segments(test_data());

    /* Successively build the tree and query intersection with the test
     * segment: */
    std::cout << "Build and check tree successively:\n";
    for (size_t i=0; i<segments.size(); ++i)
    {
        std::cout << "iteration " << i << "\n" << std::flush;

        /* Ensure that the tree is internally valid: */
        if (!tree.internally_valid())
            throw std::runtime_error("tree not valid.");

        /* Add the new segment: */
        std::cout << "   add_segment\n" << std::flush;
        tree.add_segment(segments[i]);

        /* Do the intersection test: */
        std::cout << "   intersects?\n" << std::flush;
        bool result_tree =  tree.intersects(the_segment);
        bool result_brute = brute_force_intersection(
            the_segment,
            std::vector<segment_t>(segments.begin(),
                                   segments.begin()+i+1)
        );
        if (result_tree != result_brute)
            throw std::runtime_error("Tree and brute force intersection test "
                "are not equal."
            );
    }

    if (!tree.internally_valid())
        throw std::runtime_error("Tree not internally valid.");

    /*
     * Intersection with the segment:
     */
    bool inters = tree.intersects(the_segment);
    for (auto it = tree.intersection(the_segment).begin();
         it != tree.intersection(the_segment).end(); ++it)
    {
        std::cout << "intersects at leaf " << it.leaf_id() << "\n";
    }

    /* Circle of 0.02 radius: */
    std::cout << "Intersection with r=0.02 circle.\n";
    inters = tree.intersects(the_circle);
    std::vector<size_t> I0;
    for (auto it = tree.intersection(the_circle).begin();
         it != tree.intersection(the_circle).end(); ++it)
    {
        I0.push_back(it.leaf_id());
    }
    if (I0.size() != 1 || I0.front() != 15)
        throw std::runtime_error("Segment-circle intersection at radius 0.02 "
            "is not correct."
        );

    /* Circle of 0.05 radius: */
    std::cout << "Intersection with r=0.05 circle.\n";
    I0.clear();
    for (auto it = tree.intersection(circle2).begin();
         it != tree.intersection(circle2).end(); ++it)
    {
        I0.push_back(it.leaf_id());
    }

    std::sort(I0.begin(), I0.end());
    if (I0.size() != 4 || I0[0] != 6 || I0[1] != 13 || I0[2] != 15 ||
        I0[3] != 27
    )
        throw std::runtime_error("Segment-circle intersection at radius 0.05 "
            "is not correct."
        );
}
#pragma GCC pop_options


}