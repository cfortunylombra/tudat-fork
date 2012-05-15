/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110110    K. Kumar          File created.
 *      110121    K. Kumar          Updated to comply with new protocol.
 *      110310    K. Kumar          Changed right ascension of ascending node to longitude of
 *                                  ascending node.
 *      120511    K. Kumar          Boostified unit test.
 *
 *    References
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <TudatCore/Basics/testMacros.h>

#include "Tudat/Astrodynamics/States/keplerianElements.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_keplerian_state )

//! Test if individual set- and get-functions are working correctly.
BOOST_AUTO_TEST_CASE( testKeplerianStateIndividualSetFunctions )
{
    using astrodynamics::states::KeplerianElements;
    using astrodynamics::states::semiMajorAxisIndex;
    using astrodynamics::states::eccentricityIndex;
    using astrodynamics::states::inclinationIndex;
    using astrodynamics::states::argumentOfPeriapsisIndex;
    using astrodynamics::states::longitudeOfAscendingNodeIndex;
    using astrodynamics::states::trueAnomalyIndex;
    using unit_conversions::convertDegreesToRadians;

    // Create Keplerian elements state.
    KeplerianElements keplerianElementsState;

    // Create vector of Keplerian elements: semi-major axis, eccentricity,
    // inclination, argument of periapsis, longitude of the ascending node,
    // true anomaly.
    Eigen::VectorXd keplerianElements( 6 );
    keplerianElements( semiMajorAxisIndex ) = 2.5e6;
    keplerianElements( eccentricityIndex ) = 0.1;
    keplerianElements( inclinationIndex ) = convertDegreesToRadians( 102.3 );
    keplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 125.7 );
    keplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 215.34 );
    keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );

    // Test 1: Set Keplerian elements in state object.
    keplerianElementsState.setSemiMajorAxis( keplerianElements( semiMajorAxisIndex ) );
    keplerianElementsState.setEccentricity( keplerianElements( eccentricityIndex ) );
    keplerianElementsState.setInclination( keplerianElements( inclinationIndex ) );
    keplerianElementsState.setArgumentOfPeriapsis( keplerianElements( argumentOfPeriapsisIndex ) );
    keplerianElementsState.setLongitudeOfAscendingNode( keplerianElements(
                                                       longitudeOfAscendingNodeIndex ) );
    keplerianElementsState.setTrueAnomaly( keplerianElements( trueAnomalyIndex ) );

    // Check that result from using individual get-functions gives expected results.
    BOOST_CHECK_EQUAL( keplerianElements( semiMajorAxisIndex ),
                       keplerianElementsState.getSemiMajorAxis( ) );

    BOOST_CHECK_EQUAL( keplerianElements( eccentricityIndex ),
                       keplerianElementsState.getEccentricity( ) );

    BOOST_CHECK_EQUAL( keplerianElements( inclinationIndex ),
                       keplerianElementsState.getInclination( ) );

    BOOST_CHECK_EQUAL( keplerianElements( argumentOfPeriapsisIndex ),
                       keplerianElementsState.getArgumentOfPeriapsis( ) );

    BOOST_CHECK_EQUAL( keplerianElements( longitudeOfAscendingNodeIndex ),
                       keplerianElementsState.getLongitudeOfAscendingNode( ) );

    BOOST_CHECK_EQUAL( keplerianElements( trueAnomalyIndex ),
                       keplerianElementsState.getTrueAnomaly( ) );
}

//! Test if setting state vector is working correctly.
BOOST_AUTO_TEST_CASE( testKeplerianStateSetFunction )
{
    using astrodynamics::states::KeplerianElements;
    using astrodynamics::states::semiMajorAxisIndex;
    using astrodynamics::states::eccentricityIndex;
    using astrodynamics::states::inclinationIndex;
    using astrodynamics::states::argumentOfPeriapsisIndex;
    using astrodynamics::states::longitudeOfAscendingNodeIndex;
    using astrodynamics::states::trueAnomalyIndex;
    using unit_conversions::convertDegreesToRadians;

    // Create Keplerian elements state.
    KeplerianElements keplerianElementsState;

    // Create vector of Keplerian elements: semi-major axis, eccentricity,
    // inclination, argument of periapsis, longitude of the ascending node,
    // true anomaly.
    Eigen::VectorXd keplerianElements( 6 );
    keplerianElements( semiMajorAxisIndex ) = 2.5e6;
    keplerianElements( eccentricityIndex ) = 0.1;
    keplerianElements( inclinationIndex ) = convertDegreesToRadians( 102.3 );
    keplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 125.7 );
    keplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 215.34 );
    keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );

    // Test 2: Set Keplerian elements as vector in state object.
    keplerianElementsState.state = keplerianElements;

    // Check that result from using individual get-functions gives expected results.
    TUDAT_CHECK_MATRIX_BASE( keplerianElements, keplerianElementsState.state )
            BOOST_CHECK_EQUAL( keplerianElements, keplerianElementsState.state );
}

//! Test if individual, auxilliary set- and get-functions are working correctly.
BOOST_AUTO_TEST_CASE( testKeplerianStateIndividualAuxilliarySetFunctions )
{
    using astrodynamics::states::KeplerianElements;
    using astrodynamics::states::semiMajorAxisIndex;
    using astrodynamics::states::eccentricityIndex;
    using astrodynamics::states::inclinationIndex;
    using astrodynamics::states::argumentOfPeriapsisIndex;
    using astrodynamics::states::longitudeOfAscendingNodeIndex;
    using astrodynamics::states::trueAnomalyIndex;
    using astrodynamics::states::semiLatusRectumIndex;
    using astrodynamics::states::longitudeOfPeriapsisIndex;
    using astrodynamics::states::trueLongitudeIndex;
    using astrodynamics::states::argumentOfLatitudeIndex;
    using unit_conversions::convertDegreesToRadians;

    // Create Keplerian elements state.
    KeplerianElements keplerianElementsState;

    // Create vector of Keplerian elements: semi-major axis, eccentricity,
    // inclination, argument of periapsis, longitude of the ascending node,
    // true anomaly.
    Eigen::VectorXd keplerianElements( 6 );
    keplerianElements( semiMajorAxisIndex ) = std::numeric_limits< double >::signaling_NaN( );
    keplerianElements( eccentricityIndex ) = 1.0;
    keplerianElements( inclinationIndex ) = convertDegreesToRadians( 102.3 );
    keplerianElements( argumentOfPeriapsisIndex ) = convertDegreesToRadians( 125.7 );
    keplerianElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 215.34 );
    keplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );

    // Create vector of auxilliary Keplerian elements: semi-latus rectum, longitude of periapsis,
    // true longitude, argument of latitude.
    Eigen::VectorXd auxilliaryKeplerianElements( 4 );
    auxilliaryKeplerianElements( semiLatusRectumIndex ) = 2.0e6;
    auxilliaryKeplerianElements( longitudeOfPeriapsisIndex )
            = keplerianElements( argumentOfPeriapsisIndex )
            + keplerianElements( longitudeOfAscendingNodeIndex );
    auxilliaryKeplerianElements( trueLongitudeIndex )
            = keplerianElements( argumentOfPeriapsisIndex )
            + keplerianElements( longitudeOfAscendingNodeIndex )
            + keplerianElements( trueAnomalyIndex );
    auxilliaryKeplerianElements( argumentOfLatitudeIndex )
            = keplerianElements( argumentOfPeriapsisIndex )
            + keplerianElements( trueAnomalyIndex );

    // Test 3: Set Keplerian elements in state object.
    keplerianElementsState.setSemiLatusRectum(
                auxilliaryKeplerianElements( semiLatusRectumIndex ) );
    keplerianElementsState.setEccentricity( keplerianElements( eccentricityIndex ) );
    keplerianElementsState.setInclination( keplerianElements( inclinationIndex ) );
    keplerianElementsState.setArgumentOfPeriapsis( keplerianElements( argumentOfPeriapsisIndex ) );
    keplerianElementsState.setLongitudeOfAscendingNode( keplerianElements(
                                                       longitudeOfAscendingNodeIndex ) );
    keplerianElementsState.setTrueAnomaly( keplerianElements( trueAnomalyIndex ) );

    // Check that result from using individual get-functions gives expected results.
    BOOST_CHECK_EQUAL( auxilliaryKeplerianElements( semiLatusRectumIndex ),
                       keplerianElementsState.getSemiLatusRectum( ) );

    BOOST_CHECK_EQUAL( auxilliaryKeplerianElements( longitudeOfPeriapsisIndex ),
                       keplerianElementsState.getLongitudeOfPeriapsis( ) );

    BOOST_CHECK_EQUAL( auxilliaryKeplerianElements( trueLongitudeIndex ),
                       keplerianElementsState.getTrueLongitude( ) );

    BOOST_CHECK_EQUAL( auxilliaryKeplerianElements( argumentOfLatitudeIndex ),
                       keplerianElementsState.getArgumentOfLatitude( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
