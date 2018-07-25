/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

//#define BOOST_TEST_MAIN

#include "Tudat/JsonInterface/UnitTests/unitTestSupport.h"
#include "Tudat/JsonInterface/Propagation/acceleration.h"
#include "Tudat/JsonInterface/Estimation/observation.h"
#include "Tudat/JsonInterface/Estimation/parameter.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"
#include "Tudat/SimulationSetup/EstimationSetup/createEstimatableParameters.h"
#include "Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.h"
//namespace tudat
//{

//namespace unit_tests
//{

#define INPUT( filename ) \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ).string( )

//BOOST_AUTO_TEST_SUITE( test_json_observation )

//// Test 1: sphericalHarmonicGravity
//BOOST_AUTO_TEST_CASE( test_json_acceleration_sphericalHarmonicGravity )
int main( )
{
    using namespace tudat;
    using namespace tudat::simulation_setup;
    //using namespace tudat::unit_tests;

    using namespace tudat::simulation_setup;
    using namespace tudat::observation_models;
    using namespace tudat::json_interface;
    using namespace tudat::estimatable_parameters;

    // Create ground stations from geodetic positions.
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );


    // Create list of link ends in which station is receiver and in which station is transmitter (with spacecraft other link end).
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ receiver ] = std::make_pair( "Vehicle", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = std::make_pair( "Vehicle", "" );
        stationReceiverLinkEnds.push_back( linkEnds );
    }

    std::vector< boost::shared_ptr < AccelerationSettings > > accSettings;
    accSettings.push_back( boost::make_shared< AccelerationSettings >(
                basic_astrodynamics::central_gravity ) );

    {
        nlohmann::json jsonObject;
        std::cout<<"Converting"<<std::endl;
        to_json( jsonObject, accSettings );
        std::cout<<"Converted"<<std::endl;


        std::vector< boost::shared_ptr < AccelerationSettings > > accSettingsFromJson;
        from_json( jsonObject, accSettingsFromJson );
    }

    sleep( 10000.0 );

    // Define (arbitrarily) link ends to be used for 1-way range, 1-way doppler and angular position observables
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
//    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
//    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

//    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
//    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 2 ] );

//    linkEndsPerObservable[ angular_position ].push_back( stationReceiverLinkEnds[ 2 ] );
//    linkEndsPerObservable[ angular_position ].push_back( stationTransmitterLinkEnds[ 1 ] );

    // Iterate over all observable types and associated link ends, and creatin settings for observation
    //observation_models::ObservationSettingsMap observationSettingsMap;
    std::map< LinkEnds, std::vector< boost::shared_ptr < ObservationSettings > > > observationSettingsMap;

    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;

        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            // Define bias and light-time correction settings
            boost::shared_ptr< ObservationBiasSettings > biasSettings;
            boost::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections;
//            if( currentObservable == one_way_range && i == 0)
//            {
//                biasSettings = boost::make_shared< ConstantObservationBiasSettings >(
//                            Eigen::Vector1d::Constant( 1.0 ), true );
//                std::vector< std::string > perturbingBodies = { "Mars", "Moon" };
//                lightTimeCorrections = boost::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
//                            perturbingBodies );
//            }

//            if( currentObservable == angular_position && i == 1 )
//            {
//                biasSettings = boost::make_shared< ConstantObservationBiasSettings >(
//                            Eigen::Vector2d::Constant( 1.0E-6 ), false );
//            }

//            if( currentObservable == one_way_doppler )
//            {
//                biasSettings = boost::make_shared< ConstantObservationBiasSettings >(
//                            Eigen::Vector1d::Constant( 1.0E-12 ), false );
//            }

            boost::shared_ptr< ObservationSettings > currentObservationSettings;
            if( currentObservable == one_way_doppler && i == 0 )
            {
                currentObservationSettings = boost::make_shared< OneWayDopplerObservationSettings >(
                            lightTimeCorrections,
                            boost::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                            boost::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Sun" ),
                            biasSettings );
            }
            else
            {
                currentObservationSettings = boost::make_shared< ObservationSettings >(
                            currentObservable, lightTimeCorrections, biasSettings );
            }


            // Define settings for observable, no light-time corrections, and biases for selected links
            observationSettingsMap[ currentLinkEndsList.at( i ) ].push_back( currentObservationSettings );
        }
    }


    nlohmann::json jsonObject;
    std::cout<<"Converting"<<std::endl;
    to_json( jsonObject, observationSettingsMap );
    std::cout<<"Converted"<<std::endl;

    std::string fileName = "/home/dominic/Software/numericalAstrodynamicsTudatBundle/tudatBundle/tudat/Tudat/JsonInterface/UnitTests/observationOutput.json";
    std::ofstream outputFile( fileName );
    outputFile << jsonObject.dump( 2 );
    outputFile.close( );

//    std::map< LinkEnds, std::vector< boost::shared_ptr< ObservationSettings > > > observationSettingsMapFromJson;
//    from_json( jsonObject, observationSettingsMapFromJson );

    nlohmann::json jsonObjectFromFile = parseJSONFile( fileName );

    std::map< LinkEnds, std::vector< boost::shared_ptr< ObservationSettings > > > observationSettingsMapFromFile =
            parseJSONFile< std::map< LinkEnds, std::vector< boost::shared_ptr< ObservationSettings > > > >( fileName );

    std::cout<<jsonObject<<std::endl;

    std::string linkEndString = "transmitter:(Earth, Station1); receiver:(Vehicle, )";

    LinkEnds test = boost::lexical_cast< LinkEnds >( linkEndString );
    for( auto it = test.begin( ); it != test.end( );it++ )
    {
        std::cout<<it->first<<" "<<it->second.first<<" "<<it->second.second<<std::endl;
    }

}


//BOOST_AUTO_TEST_SUITE_END( )

//} // namespace unit_tests

//} // namespace tudat
